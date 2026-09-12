# Manuscript / Supplementary consistency audit (2026-09-12)
# Scope: polished main text vs source CSVs, S0-S10, Fig 1/S0-S7, Data availability

## Status: FIXED in this packaging pass

1. **Figure S0 PRISMA was stale**
   - Was: Excluded n=30; included “+ 1 chemotherapy sensitivity cohort (n=18)”
   - Now: Excluded n=31 (ineligible design 24); included 9 nCRT (n=577) + 1 RT (n=51)
   - Script `fig_S0_prisma.R` regenerated; matches manuscript 41/31/10/628

2. **Regimen interaction Methods wording (polished tex)**
   - Was: “score × regimen interaction term in logistic response models”
   - Actual analysis (`gse87211_regimen_sensitivity.R`): `lm(z(score) ~ response * regimen)`
   - Fixed to linear model of z-scored immune score

3. **Stale pipeline artifacts after GSE213331 removal**
   - `validation_summary.csv`: re-run → config `2026-09-12-gse213331-removed-v1`, 10 cohorts only
   - `score_high_pooled.csv` / `meta_prediction_intervals.csv`: re-run → only primary_ncrt9, radiotherapy_containing10, high_coverage_ncrt7 (no chemotherapy/all11)
   - `validate_outputs.R` success message no longer mentions chemotherapy cohort

4. **Provenance tables**
   - `table_S1_response_harmonization.csv`: GSE213331 row marked REMOVED (post-treatment resections)
   - `cohort_search_record.md`: inclusion/exclusion counts unified to 41/31/10/628; GSE213331 moved to ineligible design

## Status: VERIFIED consistent (main text numbers vs source)

| Claim in polished text | Source | Match |
|---|---|---|
| Primary IFN g=0.38 (0.16–0.61), P=0.0041, I²=1%, Q=8.12, 8/9 | meta_pooled_results.csv | YES |
| Cyto g=0.26 P=0.026; NK g=0.23 P=0.054 | same | YES |
| PI IFN 0.156–0.611 vs CI 0.161–0.606 | meta_prediction_intervals.csv | YES |
| Coverage 7 cohorts IFN 0.45 / Cyto 0.32 / NK 0.26 | high_coverage_ncrt7 rows | YES |
| RT-scope IFN 0.30 P=0.030 I²=34%; Cyto P=0.058; NK P=0.106 | radiotherapy_containing10 | YES |
| GSE56699 IFN/NK/Cyto −0.32/−0.19/−0.16 | meta_all_cohorts.csv | YES |
| Top-tertile NK OR 1.98 (1.13–3.44) P=0.022 PI 0.70–5.54 | score_high_pooled.csv | YES |
| IFN OR 1.66; Cyto OR 1.57 | same | YES |
| Endpoint complete g=0.40 k=2 n=295 no PI; broader g=0.37 (0.02–0.72) P=0.043 I²=21% | meta_endpoint_*.csv | YES |
| Q_between P=0.89/0.19/0.26 | meta_endpoint_q_between.csv | YES |
| Regimen n=111/91; g=0.47/0.54; interaction P=0.70/0.79/0.19 | gse87211_regimen_*.csv | YES |
| Ayers/Reactome g=0.29/0.35/0.33 P=0.034/0.008/0.019 | ifn_robustness_pooled.csv | YES |
| STAT1/CXCL10/IFNG g=0.39/0.32/0.28 | same | YES |
| mean-z g=0.38 P=0.0066 8/9 | ifn_scoring_pooled.csv | YES |
| MCPcounter Cyto g=0.31 P=0.025; NK g=0.21 P=0.072 | mcp_package_pooled.csv | YES |
| Bootstrap ΔAUC 0.049 (−0.058–0.136); LR P=0.016 | gse87211_bootstrap_auc.csv | YES |
| AUC clin 0.600→0.655; DeLong 0.21; KRAS P=0.17; IFN|KRAS P=0.0088 | GSE87211_clinical_utility.csv | YES |
| IG3/4 vs IG1/2 Fisher P=0.012 (0.64 vs 0.24) | GSE209746_ig_groups_summary.csv | YES |
| MSS-only IG3 Fisher P=0.043 | same | YES |
| DFS 175/34; HR 1.05 (0.72–1.54) P=0.80; spline P=0.44; pCR 0 events | table_S7 / survival summary | YES |
| Sequential c-index 0.675→0.680→0.713; IFN HR 1.07 P=0.74; 1.05 P=0.81 | gse87211_sequential_survival.csv | YES |
| TCGA IFN OS P=0.70; NK PFI HR=0.70 P=0.031 | table_S7 | YES |
| LOO IFN g 0.34–0.43 P=0.0036–0.020 | meta_leave_one_out.csv | YES |
| Egger P=0.58/0.24/0.050 | meta_egger.csv | YES |
| Table 1 n and responder % | project_config / validation_summary | YES |
| Fig 1 caption numbers (0.38/0.26/0.23; Fisher 0.012; 0/32 vs 34/143) | above sources | YES |
| Data availability lists 10 GEO accessions + TCGA-READ | project cohorts | YES |
| S0–S7 and S1–S10 all defined in legends | polished tex | YES |
| Figure files present for Fig1 + S0–S7 | disk inventory | YES |

## Residual notes (not blocking; documented)

1. **Fig S1 filename**: assembled as Hallmark barplot (`hallmark_3cohorts_barplot.png`); legend says “Figure S1”. Package renames to `Figure_S1_hallmark_camera.png`.
2. **Fig S3**: `fig_S_score_high.png` → package as `Figure_S3_score_high.png`.
3. **Methods cite Table S2** for effect-size-level pooling; S2 is platform harmonization table (appropriate).
4. **AUC range 0.556–0.764** in Results uses combined_model four-cohort IFN-alone AUCs; GSE87211 clinical-subset IFN-alone AUC is 0.597 (n=197). Not contradictory (different subsets), but do not mix the two numbers.
5. **I² primary “1%”** = 1.507 rounded; abstract/table use 1%.
6. **project `manuscript.tex`** in repo was updated earlier with the same endpoint/regimen sections; polished Downloads copy is the submission master for this package (logistic wording corrected there).
7. **validate_outputs** still only validates feature CSVs, not PRISMA/manifest; PRISMA fixed by script edit this pass.

## Package contents (submission zip)

```
LARC_IFN_submission_2026-09-12/
  manuscript/
    LARC_IFN_meta_analysis_polished_2026-09-12.tex
    LARC_IFN_meta_analysis_polished_2026-09-12.pdf
    Figure_1_main.png
  figures/
    Figure_S0_prisma.png
    Figure_S1_hallmark_camera.png
    Figure_S2_diagnostics.png
    Figure_S3_score_high.png
    Figure_S4_models.png
    Figure_S5_genes15_meta.png
    Figure_S6_survival.png
    Figure_S7_ifn_robustness.png
  tables/
    Table_S1_response_harmonization.csv
    Table_S2_platform_harmonization.csv
    Table_S3_genes15_meta_pooled.csv
    Table_S3_genes15_meta_percohort.csv
    Table_S4_leave_one_out.csv
    Table_S4_egger.csv
    Table_S5_score_high_pooled.csv
    Table_S5_immune_hot_xcohort.csv
    Table_S6_confounder_adjustment.csv
    Table_S7_survival_models.csv
    Table_S8_ifn_robustness_pooled.csv
    Table_S9_ifn_scoring_mcpcounter.csv (pooled rows)
    Table_S10_endpoint_subgroups.csv
    Table_S10_endpoint_q_between.csv
    Table_S10_regimen_effects.csv
    Table_S10_regimen_interaction.csv
  provenance/
    cohort_search_record.md
    validation_summary.csv
    meta_pooled_results.csv
    CONSISTENCY_AUDIT.md
```
