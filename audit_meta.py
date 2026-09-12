import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t
from scipy import stats


COHORTS = [
    "GSE209746", "GSE35452", "GSE87211", "GSE150082", "GSE45404",
    "GSE119409", "GSE133057", "GSE53781", "GSE94104", "GSE56699",
    "GSE213331",
]
FEATURES = ["IFN", "NK", "Cyto"]


def load_registry(path="project_config.R"):
    """Derive the authoritative cohort order / analysis sets from the single R
    registry so this Python audit cannot silently drift from project_config.R."""
    import re
    from pathlib import Path

    txt = Path(path).read_text(encoding="utf-8")
    m = re.search(r"cohort\s*=\s*c\((.*?)\)", txt, flags=re.S)
    if not m:
        raise RuntimeError("cannot parse cohort list from project_config.R")
    cohorts = re.findall(r'"(GSE\d+)"', m.group(1))
    m_ncrt = re.search(r'rep\(\s*"nCRT"\s*,\s*(\d+)\s*L?\s*\)', txt)
    n_ncrt = int(m_ncrt.group(1)) if m_ncrt else len(cohorts) - 1
    ncrt = cohorts[:n_ncrt]
    m_hc = re.search(r"high_coverage_ncrt7\s*=.*?setdiff\(.*?c\((.*?)\)", txt, flags=re.S)
    hc_exclude = re.findall(r'"(GSE\d+)"', m_hc.group(1)) if m_hc else []
    analysis_sets_cfg = {
        "primary_ncrt9": ncrt,
        "chemotherapy_containing10": [c for c in cohorts if c != "GSE56699"],
        "radiotherapy_containing10": [c for c in cohorts if c != "GSE213331"],
        "all_neoadjuvant11": cohorts,
        "high_coverage_ncrt7": [c for c in ncrt if c not in hc_exclude],
    }
    return cohorts, analysis_sets_cfg


CONFIG_COHORTS, CONFIG_ANALYSIS_SETS = load_registry()
assert CONFIG_COHORTS == COHORTS, (
    f"COHORTS drift from project_config.R: audit={COHORTS} config={CONFIG_COHORTS}"
)


def hedges(values, groups):
    values = np.asarray(values, float)
    groups = np.asarray(groups)
    good = values[groups == "good"]
    poor = values[groups == "poor"]
    n_good, n_poor = len(good), len(poor)
    pooled_sd = np.sqrt(
        ((n_good - 1) * good.var(ddof=1) + (n_poor - 1) * poor.var(ddof=1))
        / (n_good + n_poor - 2)
    )
    d = (good.mean() - poor.mean()) / pooled_sd
    correction = 1 - 3 / (4 * (n_good + n_poor) - 9)
    estimate = correction * d
    se = np.sqrt(
        (n_good + n_poor) / (n_good * n_poor)
        + estimate**2 / (2 * (n_good + n_poor - 2))
    )
    return estimate, se, n_good, n_poor


def pool_random_effects(estimates, ses):
    estimates = np.asarray(estimates, float)
    ses = np.asarray(ses, float)
    variances = ses**2
    k = len(estimates)

    def restricted_nll(tau2):
        weights = 1 / (variances + tau2)
        mean = np.sum(weights * estimates) / np.sum(weights)
        return 0.5 * (
            np.sum(np.log(variances + tau2))
            + np.log(np.sum(weights))
            + np.sum(weights * (estimates - mean) ** 2)
        )

    upper = max(1.0, 100 * np.var(estimates, ddof=1), 100 * np.max(variances))
    optimum = minimize_scalar(
        restricted_nll, bounds=(0, upper), method="bounded", options={"xatol": 1e-12}
    )
    tau2 = optimum.x if restricted_nll(optimum.x) < restricted_nll(0) else 0.0
    weights = 1 / (variances + tau2)
    mean = np.sum(weights * estimates) / np.sum(weights)
    q_hk = np.sum(weights * (estimates - mean) ** 2) / (k - 1)
    scale = max(1.0, q_hk)
    se_mean = np.sqrt(scale / np.sum(weights))
    degrees_freedom = k - 1
    critical = t.ppf(0.975, degrees_freedom)
    ci_low = mean - critical * se_mean
    ci_high = mean + critical * se_mean
    p_value = 2 * t.sf(abs(mean / se_mean), degrees_freedom)

    fixed_weights = 1 / variances
    fixed_mean = np.sum(fixed_weights * estimates) / np.sum(fixed_weights)
    q = np.sum(fixed_weights * (estimates - fixed_mean) ** 2)
    i2 = max(0, (q - (k - 1)) / q) * 100 if q else 0
    prediction_critical = t.ppf(0.975, k - 2)
    prediction_half_width = prediction_critical * np.sqrt(tau2 + se_mean**2)
    return {
        "k": k,
        "est": mean,
        "lo": ci_low,
        "hi": ci_high,
        "p": p_value,
        "I2": i2,
        "tau2": tau2,
        "pi_lo": mean - prediction_half_width,
        "pi_hi": mean + prediction_half_width,
    }


frames = {}
effect_rows = []
for cohort in COHORTS:
    frame = pd.read_csv(f"features_{cohort}.csv", index_col=0)
    if cohort == "GSE87211":
        frame = frame.drop(index="GSM2325248", errors="ignore")
    frames[cohort] = frame
    for feature in FEATURES:
        estimate, se, n_good, n_poor = hedges(frame[feature], frame.response)
        effect_rows.append(
            {
                "cohort": cohort,
                "feature": feature,
                "n": len(frame),
                "n_good": n_good,
                "n_poor": n_poor,
                "g": estimate,
                "se": se,
            }
        )

effects = pd.DataFrame(effect_rows)
effects["treatment_scope"] = np.select(
    [effects.cohort.eq("GSE56699"), effects.cohort.eq("GSE213331")],
    ["neoadjuvant_radiotherapy", "neoadjuvant_chemotherapy"],
    default="nCRT",
)
effects["lo"] = effects.g - 1.96 * effects.se
effects["hi"] = effects.g + 1.96 * effects.se
analysis_sets = {
    "primary_ncrt9": [
        cohort for cohort in COHORTS if cohort not in ("GSE56699", "GSE213331")
    ],
    "chemotherapy_containing10": [cohort for cohort in COHORTS if cohort != "GSE56699"],
    "radiotherapy_containing10": [cohort for cohort in COHORTS if cohort != "GSE213331"],
    "all_neoadjuvant11": COHORTS,
    "high_coverage_ncrt7": [
        cohort
        for cohort in COHORTS
        if cohort not in ("GSE56699", "GSE213331", "GSE133057", "GSE53781")
    ],
}
assert analysis_sets == CONFIG_ANALYSIS_SETS, (
    "analysis_sets drift from project_config.R: "
    f"audit={analysis_sets} config={CONFIG_ANALYSIS_SETS}"
)

print("CONTINUOUS")
continuous_summary = []
for analysis, cohorts in analysis_sets.items():
    for feature in FEATURES:
        subset = effects[
            effects.feature.eq(feature) & effects.cohort.isin(cohorts)
        ]
        result = pool_random_effects(subset.g, subset.se)
        continuous_summary.append({
            "analysis": analysis, "feature": feature, "k": result["k"],
            "n_total": int(subset.n.sum()), "estimate": result["est"],
            "ci_lo": result["lo"], "ci_hi": result["hi"], "p": result["p"],
            "I2": result["I2"], "tau2": result["tau2"],
            "pi_lo": result["pi_lo"], "pi_hi": result["pi_hi"],
            "direction_consistent": int((subset.g > 0).sum()),
            "method": "REML + modified Hartung-Knapp",
        })
        values = " ".join(
            f"{result[key]:.8g}"
            for key in ["k", "est", "lo", "hi", "p", "I2", "tau2", "pi_lo", "pi_hi"]
        )
        print(analysis, feature, values)

print("\nGSE87211")
print(effects[effects.cohort.eq("GSE87211")].to_string(index=False))

print("\nPRIMARY NCRT DIRECTIONS, LOO, EGGER")
for feature in FEATURES:
    subset = effects[effects.feature.eq(feature) & effects.cohort.isin(analysis_sets["primary_ncrt9"])]
    print(feature, "positive", int((subset.g > 0).sum()), "of", len(subset))
    loo = []
    for omitted in subset.cohort:
        reduced = subset[~subset.cohort.eq(omitted)]
        fit = pool_random_effects(reduced.g, reduced.se)
        loo.append((omitted, fit["est"], fit["lo"], fit["hi"], fit["p"]))
    print("LOO", min(x[1] for x in loo), max(x[1] for x in loo), min(x[4] for x in loo), max(x[4] for x in loo))
    precision = 1 / subset.se.to_numpy()
    standardized = subset.g.to_numpy() / subset.se.to_numpy()
    regression = stats.linregress(precision, standardized)
    # Egger intercept test uses the intercept standard error and k-2 df.
    design = np.column_stack([np.ones(len(subset)), precision])
    beta = np.linalg.lstsq(design, standardized, rcond=None)[0]
    residual = standardized - design @ beta
    covariance = np.linalg.inv(design.T @ design) * (residual @ residual) / (len(subset) - 2)
    t_intercept = beta[0] / np.sqrt(covariance[0, 0])
    p_intercept = 2 * stats.t.sf(abs(t_intercept), df=len(subset) - 2)
    print("Egger", beta[0], t_intercept, p_intercept)

score_high_rows = []
for cohort, frame in frames.items():
    for feature in FEATURES:
        high = frame[feature] >= frame[feature].quantile(2 / 3, interpolation="linear")
        response = frame.response.eq("good")
        a = int((high & response).sum())
        b = int((high & ~response).sum())
        c = int((~high & response).sum())
        d = int((~high & ~response).sum())
        cells = np.asarray([a, b, c, d], float)
        if np.any(cells == 0):
            cells += 0.5
        log_or = np.log(cells[0] * cells[3] / (cells[1] * cells[2]))
        se = np.sqrt(np.sum(1 / cells))
        score_high_rows.append(
            {
                "cohort": cohort,
                "feature": feature,
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "logOR": log_or,
                "se": se,
            }
        )

score_high = pd.DataFrame(score_high_rows)
print("\nSCORE_HIGH")
score_high_summary = []
for analysis, cohorts in analysis_sets.items():
    for feature in FEATURES:
        subset = score_high[
            score_high.feature.eq(feature) & score_high.cohort.isin(cohorts)
        ]
        result = pool_random_effects(subset.logOR, subset.se)
        score_high_summary.append({
            "analysis": analysis, "feature": feature, "k": result["k"],
            "n_total": int(sum(subset[["a", "b", "c", "d"]].sum(axis=1))),
            "OR": np.exp(result["est"]), "ci_lo": np.exp(result["lo"]),
            "ci_hi": np.exp(result["hi"]), "p": result["p"], "I2": result["I2"],
            "tau2_log": result["tau2"], "pi_lo": np.exp(result["pi_lo"]),
            "pi_hi": np.exp(result["pi_hi"]), "method": "REML + modified Hartung-Knapp",
        })
        transformed = {
            key: np.exp(result[key]) if key in {"est", "lo", "hi", "pi_lo", "pi_hi"} else result[key]
            for key in result
        }
        values = " ".join(
            f"{transformed[key]:.8g}"
            for key in ["k", "est", "lo", "hi", "p", "I2", "tau2", "pi_lo", "pi_hi"]
        )
        print(analysis, feature, values)

print("\nGSE87211 score-high")
print(score_high[score_high.cohort.eq("GSE87211")].to_string(index=False))

print("\nGENES15 (existing per-cohort effects; method audit only)")
gene_rows = pd.read_csv("genes15_meta_percohort.csv")
gene_summary = []
for gene, subset in gene_rows.dropna(subset=["g", "se"]).groupby("gene"):
    if len(subset) < 3:
        continue
    result = pool_random_effects(subset.g, subset.se)
    gene_summary.append({"gene": gene, **result})
gene_summary = pd.DataFrame(gene_summary)
# Benjamini-Hochberg adjustment using sorted indices.  Keep this as the single
# implementation so the audit cannot transiently report an incorrect FDR.
order = np.argsort(gene_summary.p.to_numpy())
ranked = gene_summary.p.to_numpy()[order] * len(gene_summary) / np.arange(1, len(gene_summary) + 1)
adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
gene_summary["fdr"] = np.nan
gene_summary.loc[gene_summary.index[order], "fdr"] = np.minimum(adjusted, 1)
print(gene_summary.sort_values("p")[["gene", "k", "est", "lo", "hi", "p", "fdr", "I2", "pi_lo", "pi_hi"]].to_string(index=False))

print("\nCSV_META_EFFECTS")
print(effects[["cohort", "treatment_scope", "feature", "n", "n_good", "n_poor", "g", "se", "lo", "hi"]].to_csv(index=False), end="")
print("CSV_META_POOLED")
print(pd.DataFrame(continuous_summary).to_csv(index=False), end="")
print("CSV_SCORE_HIGH_POOLED")
print(pd.DataFrame(score_high_summary).to_csv(index=False), end="")

# Verify the checked-in audit snapshot against the independent calculations above.
saved_effects = pd.read_csv("meta_all_cohorts.csv")
effect_keys = ["cohort", "treatment_scope", "feature", "n", "n_good", "n_poor"]
assert saved_effects[effect_keys].equals(
    effects[effect_keys].reset_index(drop=True)
)
for column in ["g", "se", "lo", "hi"]:
    assert np.allclose(saved_effects[column], effects[column], rtol=1e-12, atol=1e-12)

saved_pooled = pd.read_csv("meta_pooled_results.csv")
calculated_pooled = pd.DataFrame(continuous_summary)
assert saved_pooled[["analysis", "feature", "k", "n_total", "direction_consistent"]].equals(
    calculated_pooled[["analysis", "feature", "k", "n_total", "direction_consistent"]]
)
# estimate/CI/p/tau2/PI depend on the REML tau2 optimizer; R optimize() and
# scipy minimize_scalar legitimately differ (measured max abs diff: estimate
# 3e-5, CI 7e-5, p 3e-5, tau2 3e-5; PI bounds ~3.3e-4 because sqrt(tau2+se^2)
# amplifies the tau2 gap). atol=5e-4 covers that engine gap while remaining
# ~20x tighter than the 2-decimal reporting precision, so any real divergence
# still fails.
for column in ["estimate", "ci_lo", "ci_hi", "p", "tau2", "pi_lo", "pi_hi"]:
    assert np.allclose(saved_pooled[column], calculated_pooled[column],
                       rtol=1e-6, atol=5e-4), column
# I2 is closed-form from fixed-effect weights and must agree to near-machine precision.
assert np.allclose(saved_pooled["I2"], calculated_pooled["I2"], rtol=1e-9, atol=1e-9)

saved_high = pd.read_csv("score_high_pooled.csv")
calculated_high = pd.DataFrame(score_high_summary)
assert saved_high[["analysis", "feature", "k", "n_total"]].equals(
    calculated_high[["analysis", "feature", "k", "n_total"]]
)
# OR/CI/p/tau2/PI inherit optimizer-level differences via exp() of the pooled
# logOR (measured max abs diff <=1.3e-4); same engine-gap tolerance as above.
for column in ["OR", "ci_lo", "ci_hi", "p", "tau2_log", "pi_lo", "pi_hi"]:
    assert np.allclose(saved_high[column], calculated_high[column],
                       rtol=1e-6, atol=5e-4), column
assert np.allclose(saved_high["I2"], calculated_high["I2"], rtol=1e-9, atol=1e-9)

print("AUDIT_ASSERTIONS_PASS")
