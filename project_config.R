# Shared project registry and schema for the LARC transcriptomic meta-analysis.
# Keep cohort counts and analysis-set membership in one place so that validation,
# meta-analysis and the orchestration script cannot silently drift apart.
#
# 2026-09-12: GSE213331 REMOVED from all analysis sets. Its 18 response-labelled
# samples are POST-nCRT resected tissues (9 pCR + 9 npCR, per GEO
# !Sample_description and !Series_overall_design); the only pretreatment
# samples (9 biopsies) carry no response label. Using resected tissue as a
# pretreatment response-biomarker cohort is a design violation. The features
# file is retained on disk for provenance but is excluded from every analysis.

LARC_CONFIG <- local({
  cohorts <- data.frame(
    cohort = c(
      "GSE209746", "GSE35452", "GSE87211", "GSE150082", "GSE45404",
      "GSE119409", "GSE133057", "GSE53781", "GSE94104", "GSE56699"
    ),
    treatment_scope = c(
      rep("nCRT", 9L), "neoadjuvant_radiotherapy"
    ),
    n = c(93L, 46L, 202L, 39L, 42L, 56L, 33L, 26L, 40L, 51L),
    n_good = c(26L, 24L, 35L, 16L, 19L, 15L, 13L, 10L, 11L, 23L),
    n_poor = c(67L, 22L, 167L, 23L, 23L, 41L, 20L, 16L, 29L, 28L),
    stringsAsFactors = FALSE
  )

  # These are the columns required in every cohort-level feature table.
  # The first CSV column is the sample identifier written by write.csv().
  # IGF2/L1CAM are retained for gene-level or model analyses, but incomplete
  # historical platform coverage means they are not required to be finite.
  required_feature_columns <- c("IFN", "IGF2", "L1CAM", "NK", "Cyto", "response")
  meta_feature_columns <- c("IFN", "NK", "Cyto")
  optional_gene_columns <- c("IGF2", "L1CAM")

  analysis_sets <- list(
    primary_ncrt9 = cohorts$cohort[cohorts$treatment_scope == "nCRT"],
    radiotherapy_containing10 = cohorts$cohort[
      cohorts$treatment_scope %in% c("nCRT", "neoadjuvant_radiotherapy")
    ],
    high_coverage_ncrt7 = setdiff(
      cohorts$cohort[cohorts$treatment_scope == "nCRT"],
      c("GSE133057", "GSE53781")
    )
  )

  list(
    version = "2026-09-12-gse213331-removed-v1",
    cohorts = cohorts,
    required_feature_columns = required_feature_columns,
    meta_feature_columns = meta_feature_columns,
    optional_gene_columns = optional_gene_columns,
    analysis_sets = analysis_sets
  )
})

stopifnot(
  is.data.frame(LARC_CONFIG$cohorts),
  identical(
    LARC_CONFIG$cohorts$n,
    LARC_CONFIG$cohorts$n_good + LARC_CONFIG$cohorts$n_poor
  ),
  sum(LARC_CONFIG$cohorts$n) == 628L,
  sum(LARC_CONFIG$cohorts$n[LARC_CONFIG$cohorts$treatment_scope == "nCRT"]) == 577L,
  sum(LARC_CONFIG$cohorts$n[
    LARC_CONFIG$cohorts$treatment_scope == "neoadjuvant_radiotherapy"
  ]) == 51L
)
