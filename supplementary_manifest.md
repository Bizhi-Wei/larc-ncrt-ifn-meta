# Supplementary output manifest

This manifest maps the manuscript legends (2026-09-11 GSE56699 audit revision) to their source scripts and files. Final generated status is based on the successful 37-step full run recorded in `run_all_20260911_121136.log`.

| Item | Intended content | Source script | Output files | Status |
|---|---|---|---|---|
| Figure S0 | PRISMA-style cohort flow | `fig_S0_prisma.R` | `fig_S0_prisma.png`, `fig_S0_prisma.svg` | PASS: regenerated 2026-09-11 |
| Figure S1 | Hallmark pathway consistency across four initial cohorts | `hallmark_3cohorts.R` (legacy name retained) | `hallmark_3cohorts_barplot.png`, `.pdf`, `.svg`; `hallmark_4cohorts_camera.csv` | PASS: regenerated and visually checked 2026-09-11 |
| Figure S2 | Leave-one-cohort-out, centred funnel plots, Baujat influence | `meta_diagnostics.R`, `meta_baujat.R`, `fig_S2_diagnostics.R` | `fig_S2_diagnostics.png`, `.pdf`, `.svg` | PASS: regenerated and visually checked 2026-09-11 |
| Figure S3 | Top-tertile score-high sensitivity analysis | `immune_hot_xcohort.R` | `fig_S_score_high.png` | PASS: regenerated 2026-09-11 |
| Figure S4 | IFN ROC and fixed-coefficient transfer model | `IFN_signature_ROC_MSI.R`, `combined_model.R`, `fig_supplement_assemblies.R` | `fig_S4_models.png`, `.pdf`, `.svg` | PASS: regenerated and visually checked 2026-09-11 |
| Figure S5 | Fifteen-gene REML-mKH meta-analysis | `genes15_meta.R` | `fig_S_genes15_meta.png` | PASS: regenerated 2026-09-11; no BH-significant gene |
| Figure S6 | GSE87211 spline and TCGA-READ survival context | `GSE87211_survival_v2.R`, `TCGA_READ_survival.R`, `fig_supplement_assemblies.R` | `fig_S6_survival.png`, `.pdf`, `.svg` | PASS: regenerated and visually checked 2026-09-11 |
| Figure S7 | IFN-signature robustness across independent gene sets | `ifn_robustness.R` | `fig_S_ifn_robustness.png` | PASS: regenerated 2026-09-11 |
| Table S1 | Response-definition harmonization per cohort | manual curation, verified against analysis_log.md | `table_S1_response_harmonization.csv` | Updated for GSE56699 patient-level mapping |
| Table S2 | Platform and preprocessing harmonization chain | manual curation from per-cohort scripts | `table_S2_platform_harmonization.csv` | Updated for official GSE56699 inputs |
| Table S3 | Fifteen-gene per-cohort and pooled effects | `genes15_meta.R` | `genes15_meta_percohort.csv`, `genes15_meta_pooled.csv` | PASS: regenerated 2026-09-11 |
| Table S4 | Leave-one-cohort-out results | `meta_diagnostics.R` | `meta_leave_one_out.csv`, `meta_egger.csv` | PASS: regenerated 2026-09-11 |
| Table S5 | Per-cohort top-tertile tables and pooled ORs | `immune_hot_xcohort.R` | `immune_hot_xcohort.csv`, `score_high_pooled.csv` | PASS: regenerated 2026-09-11 |
| Table S6 | MSI, stromal, WES and KRAS covariate analyses | Component analysis scripts plus `assemble_supplementary_tables.R` | `table_S6_confounder_adjustment.csv`; `outputs/01a08bc2-2052-7d60-8183-6fe862fb032d/supplementary_tables_S6_S7.xlsx`, sheet `Table S6` | PASS: 30 rows across 4 sections; regenerated 2026-09-11 |
| Table S7 | GSE87211 and TCGA-READ survival models | Survival scripts plus `assemble_supplementary_tables.R` | `table_S7_survival_models.csv`; same workbook, sheet `Table S7` | PASS: 17 rows across 2 cohorts and 3 endpoints; regenerated 2026-09-11 |
| Table S8 | IFN-signature robustness: alternative signatures and core genes | `ifn_robustness.R` | `ifn_robustness_percohort.csv`, `ifn_robustness_pooled.csv`, `ifn_robustness_correlations.csv` | PASS: regenerated 2026-09-11 |
| Table S9 | Implementation sensitivity (mean-z scoring; official MCPcounter) | `ifn_scoring_sensitivity.R`, `mcp_package_sensitivity.R` | `ifn_scoring_pooled.csv`, `ifn_scoring_sensitivity.csv`, `mcp_package_pooled.csv`, `mcp_package_sensitivity.csv` | PASS: 2026-09-12 |
| Table S10 | Response-definition subgroups and GSE87211 regimen sensitivity | `meta_endpoint_subgroups.R`, `gse87211_regimen_sensitivity.R` | `meta_endpoint_subgroups.csv`, `meta_endpoint_q_between.csv`, `meta_endpoint_effects_classified.csv`, `gse87211_regimen_counts.csv`, `gse87211_regimen_effects.csv`, `gse87211_regimen_interaction.csv` | PASS: 2026-09-12; complete-response subgroup k=2, no PI by design |

Convention note (uniform across all figures/tables since 2026-09-09): positive Hedges g = higher score in responders.

The systematic-search record is `cohort_search_record.md`. The updated pipeline, including Table S6/S7 assembly, completed successfully and was cross-checked against `validation_summary.csv`, `audit_meta.py`, and `run_all_20260911_121136.log`.
