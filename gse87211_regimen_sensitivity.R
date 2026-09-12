# ==============================
# GSE87211 regimen sensitivity (within-cohort)
# Endpoint: pCR = ypT0N0 (unchanged from primary coding)
# Strata: fluoropyrimidine-CRT (5-FU + RT)
#         vs oxaliplatin-containing CRT (5-FU + Oxa ± Cetuximab + RT)
# Cetuximab arm (n small) is merged into the oxaliplatin-containing stratum.
# Outputs:
#   gse87211_regimen_counts.csv
#   gse87211_regimen_effects.csv   (stratified Hedges g)
#   gse87211_regimen_interaction.csv
# ==============================

suppressMessages({
  library(GEOquery)
})
source("meta_utils.R", encoding = "UTF-8")

features <- read.csv("features_GSE87211.csv", stringsAsFactors = FALSE)
stopifnot(nrow(features) == 202L)
rownames(features) <- features[[1]]
features$sample_id <- rownames(features)

eset <- GEOquery:::parseGSEMatrix(
  "GSE87211_series_matrix.txt.gz",
  AnnotGPL = FALSE, getGPL = FALSE
)$eset
clin <- pData(eset)

therapy_col <- grep(
  "preoperative radiochemotherapy",
  colnames(clin), value = TRUE, ignore.case = TRUE
)[1]
stopifnot(!is.na(therapy_col))
tissue_col <- grep("^tissue", colnames(clin), value = TRUE)[1]
stopifnot(sum(clin[[tissue_col]] == "rectal tumor") == 203L)

therapy_raw <- as.character(clin[[therapy_col]])
regimen3 <- ifelse(
  grepl("Cetuximab", therapy_raw, ignore.case = TRUE), "oxa_cetux",
  ifelse(grepl("Oxaliplatin", therapy_raw, ignore.case = TRUE), "oxa",
         ifelse(grepl("5-FU", therapy_raw, ignore.case = TRUE), "fluoro", NA_character_))
)
stopifnot(all(!is.na(regimen3)))

regimen_bin <- ifelse(regimen3 == "fluoro", "fluoropyrimidine_CRT",
                      "oxaliplatin_containing_CRT")

# pCR from the same three-state rule as GSE87211_analysis.R
ypT <- suppressWarnings(as.numeric(as.character(
  clin[[grep("depth of invasion after", colnames(clin), value = TRUE)[1]]]
)))
ypN <- suppressWarnings(as.numeric(as.character(
  clin[[grep("lymph node metastasis after", colnames(clin), value = TRUE)[1]]]
)))
pcr <- ifelse(ypT == 0 & ypN == 0, "good",
              ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0),
                     "poor", NA_character_))

tum_ids <- colnames(eset)[clin[[tissue_col]] == "rectal tumor"]
stopifnot(all(rownames(features) %in% tum_ids))

meta <- data.frame(
  sample_id = tum_ids,
  regimen3 = regimen3[clin[[tissue_col]] == "rectal tumor"],
  regimen = regimen_bin[clin[[tissue_col]] == "rectal tumor"],
  response = pcr[clin[[tissue_col]] == "rectal tumor"],
  therapy_raw = therapy_raw[clin[[tissue_col]] == "rectal tumor"],
  stringsAsFactors = FALSE
)
meta <- merge(meta, features[, c("sample_id", "IFN", "NK", "Cyto")],
              by = "sample_id", all.x = TRUE)
meta <- meta[meta$sample_id %in% rownames(features), ]
stopifnot(nrow(meta) == 202L)
stopifnot(identical(sort(meta$response), sort(features$response)))

# Primary-effect samples only (pCR known)
ok <- meta$response %in% c("good", "poor")
dat <- meta[ok, ]
stopifnot(nrow(dat) == 202L)

counts <- as.data.frame.matrix(table(dat$regimen, dat$response))
counts$regimen <- rownames(counts)
counts$n <- counts$good + counts$poor
counts$pcr_rate <- counts$good / counts$n
counts$regimen3_levels <- vapply(counts$regimen, function(r) {
  paste(sort(unique(dat$regimen3[dat$regimen == r])), collapse = " | ")
}, character(1))
write.csv(counts, "gse87211_regimen_counts.csv", row.names = FALSE)

features3 <- c("IFN", "NK", "Cyto")
effect_rows <- list()
for (feature in features3) {
  for (reg in unique(dat$regimen)) {
    sub <- dat[dat$regimen == reg, ]
    eff <- hedges_g(sub[[feature]], sub$response)
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      feature = feature,
      regimen = reg,
      n = nrow(sub),
      n_good = unname(eff["n_good"]),
      n_poor = unname(eff["n_poor"]),
      g = unname(eff["g"]),
      se = unname(eff["se"]),
      lo = unname(eff["g"] - 1.96 * eff["se"]),
      hi = unname(eff["g"] + 1.96 * eff["se"])
    )
  }
}
effects <- do.call(rbind, effect_rows)
write.csv(effects, "gse87211_regimen_effects.csv", row.names = FALSE)

# Effect modification: score ~ response + regimen + response:regimen
# Interaction coefficient tests whether the good-vs-poor mean difference
# differs by regimen stratum (equal-variance linear model on z-scored feature).
interaction_rows <- list()
for (feature in features3) {
  d <- dat
  d$score <- as.numeric(scale(d[[feature]]))
  d$response <- factor(d$response, levels = c("poor", "good"))
  d$regimen <- factor(d$regimen, levels = c("fluoropyrimidine_CRT",
                                            "oxaliplatin_containing_CRT"))
  fit <- lm(score ~ response * regimen, data = d)
  sm <- summary(fit)$coefficients
  # interaction term name: responsegood:regimenoxaliplatin_containing_CRT
  inter_name <- grep("responsegood:regimen", rownames(sm), value = TRUE)
  stopifnot(length(inter_name) == 1L)
  main_name <- "responsegood"
  interaction_rows[[length(interaction_rows) + 1L]] <- data.frame(
    feature = feature,
    n = nrow(d),
    beta_good_vs_poor = sm[main_name, "Estimate"],
    se_good_vs_poor = sm[main_name, "Std. Error"],
    p_good_vs_poor = sm[main_name, "Pr(>|t|)"],
    beta_interaction = sm[inter_name, "Estimate"],
    se_interaction = sm[inter_name, "Std. Error"],
    p_interaction = sm[inter_name, "Pr(>|t|)"],
    model = "z(score) ~ response * regimen (fluoropyrimidine_CRT reference)"
  )
}
interaction <- do.call(rbind, interaction_rows)
write.csv(interaction, "gse87211_regimen_interaction.csv", row.names = FALSE)

cat("\n===== GSE87211 regimen counts =====\n")
print(counts, row.names = FALSE)
cat("\n===== Stratified Hedges g =====\n")
print(effects, row.names = FALSE, digits = 3)
cat("\n===== Response x regimen interaction =====\n")
print(interaction, row.names = FALSE, digits = 3)
cat("\nNote: cetuximab-containing samples merged into oxaliplatin-containing CRT.\n")
