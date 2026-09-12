# LARC nCRT interferon meta-analysis (code + derived tables)

Reproducible analysis pipeline for:

> **Pretreatment interferon activation and cytotoxic immune enrichment are associated with response to neoadjuvant chemoradiotherapy in locally advanced rectal cancer: a systematic multi-cohort transcriptomic meta-analysis**

Author: Bizhi Wei (Puai Medical College, Shaoyang University)  
ORCID: [0009-0008-9481-3024](https://orcid.org/0009-0008-9481-3024)

## Scope

| Set | Composition | n |
|---|---|---:|
| **Primary** | 9 pretreatment nCRT cohorts | 577 |
| Treatment-scope sensitivity | Primary + GSE56699 (preoperative radiotherapy only) | 628 |
| High-coverage sensitivity | Primary minus GSE133057 and GSE53781 (<80% probe annotation) | 518 |
| Prognostic context only | TCGA-READ (not nCRT) | 164 |

Primary result (REML + modified Hartung–Knapp; positive Hedges *g* = higher score in responders):

| Feature | *g* (95% CI) | *P* | *I*² | Direction |
|---|---|---:|---:|---|
| IFN (primary) | 0.38 (0.16–0.61) | 0.0041 | 2% | 8/9 |
| Cytotoxic lymphocytes | 0.26 (0.04–0.49) | 0.026 | 0% | 7/9 |
| NK cells | 0.23 (−0.005–0.47) | 0.054 | 0% | 7/9 |

GSE213331 was assessed and **excluded** (response-labelled samples are post-treatment resections). It is not part of any analysis set.

## Repository layout

| Path | Role |
|---|---|
| `run_all.R` | End-to-end pipeline entry (preflight + staged scripts) |
| `project_config.R` | Single source of truth for cohorts, *n*, analysis sets |
| `meta_utils.R` | Hedges *g*, REML–modified Hartung–Knapp, DerSimonian–Laird |
| `validate_outputs.R` | Pre-meta assertions on feature tables |
| `meta_*.R`, `*_sensitivity.R` | Primary meta, diagnostics, endpoint/regimen/implementation sensitivities |
| `GSE*.R`, `new_cohorts_process*.R` | Per-cohort processing and feature construction |
| `figures_*.R`, `fig_S*.R` | Main and supplementary figures |
| `features_GSE*.csv` | Within-cohort IFN/NK/Cyto (+ optional genes) and response labels |
| `meta_pooled_results.csv` and related CSVs | Derived tables used in the manuscript |
| `cohort_search_record.md` | PRISMA-style search and exclusion record |
| `analysis_log.md` | Dated decision log |
| `IFN_model_card.md`, `model_inputs_sha256.txt` | Locked scoring definition and input hashes |
| `manuscript/` | Submission LaTeX source + Figure 1 |

Large public GEO/ArrayExpress matrices and serialized R objects are **not** stored here; download them from the accessions listed below and place them as documented in `run_all.R` preflight inputs.

## Data sources (public)

- GEO: GSE209746, GSE35452, GSE87211, GSE150082, GSE45404, GSE119409, GSE133057, GSE53781, GSE94104, GSE56699
- GSE56699 processed matrix/SDRF/ADF: ArrayExpress/BioStudies **E-GEOD-56699**
- TCGA-READ: UCSC Xena (star counts, phenotype, survival)

## Environment

- R 4.6.1 (Windows used in development; Linux/macOS should work)
- Bioconductor/packages: GEOquery, GSVA, Biobase, DESeq2, limma, pROC, survival, ggplot2, gridExtra, org.Hs.eg.db, illuminaHumanv3.db, scales, svglite, png
- Optional audit: Python 3 + numpy/scipy/pandas (`audit_meta.py`)

```bash
Rscript run_all.R --preflight-only   # dependency/input check
Rscript run_all.R                    # full pipeline
python audit_meta.py                 # independent numeric audit (after CSVs exist)
```

## Reproducibility notes

- Cohort membership, sample counts and analysis sets are asserted in `project_config.R` (version string recorded in `validation_summary.csv`).
- Immune scores are ssGSEA of published MCP-counter **marker sets** (not the MCPcounter package abundance estimator, except in the dedicated package sensitivity analysis).
- Positive Hedges *g* always means higher score in responders/poor responders coding is harmonized in `features_GSE*.csv`.
- The complete-response endpoint subgroup has *k*=2: pooled estimate is exploratory and a prediction interval is not reported (*t* df = *k*−2 = 0).

## Citation

Please cite the manuscript and this repository (see `CITATION.cff`).  
After the first Zenodo–GitHub archival release, replace the repository URL in the manuscript Data availability statement with the Zenodo DOI.

## License

Code: MIT (see `LICENSE`).  
Derived numeric tables in this repository: released under the same terms for reuse with attribution.
