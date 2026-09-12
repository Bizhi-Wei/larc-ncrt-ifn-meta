# Cohort-level effect sizes and random-effects meta-analysis.
# Primary analysis: nine nCRT cohorts.  GSE56699 (radiotherapy-only) and
# GSE213331 removed 2026-09-12 (post-treatment tissue; no eligible labels).

suppressMessages({
  library(ggplot2)
  library(gridExtra)
})
source("project_config.R", encoding = "UTF-8")
source("meta_utils.R")

cohort_label <- c(
  GSE209746 = "GSE209746 (RNA-seq, CR/iCR)",
  GSE35452  = "GSE35452 (Affymetrix, TRG)",
  GSE87211  = "GSE87211 (Agilent, ypT0N0)",
  GSE150082 = "GSE150082 (Agilent, TRG)",
  GSE45404  = "GSE45404 (Affymetrix, Mandard)",
  GSE119409 = "GSE119409 (Affymetrix, sens/resist)",
  GSE133057 = "GSE133057 (Illumina, AJCC 0-1/2-3)",
  GSE53781  = "GSE53781 (CodeLink, TRG1-2/3-5)",
  GSE94104  = "GSE94104 (Illumina, TRG3/1-2)",
  GSE56699  = "GSE56699 (Illumina DASL, radiotherapy, Mandard 1-2/3-5)"
)
cohort_scope <- setNames(LARC_CONFIG$cohorts$treatment_scope,
                         LARC_CONFIG$cohorts$cohort)

files <- list.files(pattern = "^features_GSE[0-9]+\\.csv$")
available_cohorts <- sub("features_(GSE[0-9]+)\\.csv", "\\1", files)
cohorts <- LARC_CONFIG$cohorts$cohort[LARC_CONFIG$cohorts$cohort %in% available_cohorts]
stopifnot(identical(cohorts, LARC_CONFIG$cohorts$cohort))
cat("Available cohorts:", paste(cohorts, collapse = ", "), "\n")

features <- c("IFN", "NK", "Cyto")
rows <- list()
for (cohort in cohorts) {
  dat <- read.csv(paste0("features_", cohort, ".csv"), stringsAsFactors = FALSE)
  stopifnot(all(c("response", features) %in% colnames(dat)))
  stopifnot(all(dat$response %in% c("good", "poor")))
  expected_row <- LARC_CONFIG$cohorts[LARC_CONFIG$cohorts$cohort == cohort, , drop = FALSE]
  stopifnot(nrow(dat) == expected_row$n,
            sum(dat$response == "good") == expected_row$n_good,
            sum(dat$response == "poor") == expected_row$n_poor)
  stopifnot(all(vapply(dat[features], function(x) all(is.finite(x)), logical(1))))
  for (feature in features) {
    effect <- hedges_g(dat[[feature]], dat$response)
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort,
      treatment_scope = unname(cohort_scope[cohort]),
      feature = feature,
      n = nrow(dat),
      n_good = unname(effect["n_good"]),
      n_poor = unname(effect["n_poor"]),
      g = unname(effect["g"]),
      se = unname(effect["se"])
    )
  }
}
effects <- do.call(rbind, rows)
effects$lo <- effects$g - 1.96 * effects$se
effects$hi <- effects$g + 1.96 * effects$se
write.csv(effects, "meta_all_cohorts.csv", row.names = FALSE)

analysis_sets <- LARC_CONFIG$analysis_sets

pooled_rows <- list()
dl_rows <- list()
for (set_name in names(analysis_sets)) {
  for (feature in features) {
    dat <- effects[
      effects$cohort %in% analysis_sets[[set_name]] & effects$feature == feature,
    ]
    fit <- pool_random_effects(dat$g, dat$se)
    dl_fit <- pool_dl(dat$g, dat$se)
    pooled_rows[[length(pooled_rows) + 1L]] <- data.frame(
      analysis = set_name, feature = feature, k = fit$k, n_total = sum(dat$n),
      estimate = fit$g, ci_lo = fit$lo, ci_hi = fit$hi, p = fit$p,
      I2 = fit$I2, tau2 = fit$tau2, Q = fit$Q, Qp = fit$Qp,
      pi_lo = fit$pi_lo, pi_hi = fit$pi_hi,
      direction_consistent = sum(dat$g > 0), method = fit$method
    )
    dl_rows[[length(dl_rows) + 1L]] <- data.frame(
      analysis = set_name, feature = feature, k = dl_fit$k, n_total = sum(dat$n),
      estimate = dl_fit$g, ci_lo = dl_fit$lo, ci_hi = dl_fit$hi, p = dl_fit$p,
      I2 = dl_fit$I2, tau2 = dl_fit$tau2,
      direction_consistent = sum(dat$g > 0), method = dl_fit$method
    )
  }
}
pooled <- do.call(rbind, pooled_rows)
pooled_dl <- do.call(rbind, dl_rows)
write.csv(pooled, "meta_pooled_results.csv", row.names = FALSE)
write.csv(pooled_dl, "meta_pooled_dl_sensitivity.csv", row.names = FALSE)

cat("\n===== Primary REML + modified Hartung-Knapp results =====\n")
print(pooled[pooled$analysis == "primary_ncrt9", ], row.names = FALSE, digits = 3)

make_forest <- function(dat, title) {
  fit <- pool_random_effects(dat$g, dat$se)
  n <- nrow(dat)
  dat$y <- n - seq_len(n) + 1.8
  x_min <- min(dat$lo, fit$pi_lo, -0.6) - 0.1
  x_max <- max(dat$hi, fit$pi_hi, 0.6) + 0.15
  diamond <- data.frame(
    x = c(fit$lo, fit$g, fit$hi, fit$g), y = c(1, 1.3, 1, 0.7)
  )
  labels <- paste0(cohort_label[dat$cohort], " n=", dat$n)
  ggplot(dat) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0.4, ymax = 1.5,
             fill = "grey95", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi, y = y), orientation = "y",
                  width = 0.18, color = "grey40") +
    geom_point(aes(x = g, y = y), shape = 15, size = 2.2, color = "#0072B2") +
    geom_polygon(data = diamond, aes(x = x, y = y), fill = "#D55E00") +
    scale_y_continuous(breaks = c(dat$y, 1), labels = c(labels, "Pooled (REML-mKH)")) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    labs(
      title = sprintf("%s: g=%.2f (%.2f-%.2f)\nP=%.3g, I2=%.0f%%",
                      title, fit$g, fit$lo, fit$hi, fit$p, fit$I2),
      x = "Hedges g (responders vs non-responders)", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9.5, face = "bold"))
}

primary <- effects[effects$cohort %in% analysis_sets$primary_ncrt9, ]
plots <- lapply(features, function(feature) {
  make_forest(primary[primary$feature == feature, ], feature)
})
forest_grob <- arrangeGrob(grobs = plots, ncol = 3)
png("fig_S_forest_ncrt_cohorts.png", width = 4200, height = 2200, res = 300)
grid::grid.draw(forest_grob)
dev.off()
ggsave("fig_S_forest_ncrt_cohorts.pdf", forest_grob, width = 14, height = 7.3,
       device = grDevices::cairo_pdf)

cat("\nOutputs: meta_all_cohorts.csv / meta_pooled_results.csv / ",
    "meta_pooled_dl_sensitivity.csv / fig_S_forest_ncrt_cohorts.png/.pdf\n", sep = "")
cat("===== Complete =====\n")
