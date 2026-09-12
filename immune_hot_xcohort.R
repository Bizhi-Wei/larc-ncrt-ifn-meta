# Cohort-specific top-tertile score analysis.
# This is a dichotomized sensitivity analysis, not an independent biomarker validation.

suppressMessages({
  library(ggplot2)
  library(gridExtra)
})
source("project_config.R", encoding = "UTF-8")
source("meta_utils.R")

cohorts <- LARC_CONFIG$cohorts$cohort
cohort_scope <- setNames(LARC_CONFIG$cohorts$treatment_scope, cohorts)
features <- c("IFN", "NK", "Cyto")

rows <- list()
for (cohort in cohorts) {
  dat <- read.csv(paste0("features_", cohort, ".csv"), stringsAsFactors = FALSE)
  for (feature in features) {
    score <- dat[[feature]]
    score_high <- score >= stats::quantile(score, 2 / 3, na.rm = TRUE)
    responder <- dat$response == "good"
    a <- sum(score_high & responder)
    b <- sum(score_high & !responder)
    c_cell <- sum(!score_high & responder)
    d <- sum(!score_high & !responder)
    table_cells <- c(a, b, c_cell, d)
    corrected <- if (any(table_cells == 0)) table_cells + 0.5 else table_cells
    log_or <- log((corrected[1] * corrected[4]) / (corrected[2] * corrected[3]))
    se <- sqrt(sum(1 / corrected))
    fisher <- stats::fisher.test(matrix(table_cells, nrow = 2, byrow = TRUE))
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort, treatment_scope = unname(cohort_scope[cohort]), feature = feature,
      n_high = a + b, responders_high = a, response_rate_high = a / (a + b),
      n_rest = c_cell + d, responders_rest = c_cell,
      response_rate_rest = c_cell / (c_cell + d),
      OR_fisher = unname(fisher$estimate), fisher_p = fisher$p.value,
      logOR = log_or, se = se
    )
  }
}
results <- do.call(rbind, rows)
write.csv(results, "immune_hot_xcohort.csv", row.names = FALSE)

analysis_sets <- LARC_CONFIG$analysis_sets

pooled_rows <- list()
for (set_name in names(analysis_sets)) {
  for (feature in features) {
    dat <- results[
      results$cohort %in% analysis_sets[[set_name]] & results$feature == feature,
    ]
    fit <- pool_random_effects(dat$logOR, dat$se)
    pooled_rows[[length(pooled_rows) + 1L]] <- data.frame(
      analysis = set_name, feature = feature, k = fit$k,
      n_total = sum(dat$n_high + dat$n_rest),
      OR = exp(fit$g), ci_lo = exp(fit$lo), ci_hi = exp(fit$hi), p = fit$p,
      I2 = fit$I2, tau2_log = fit$tau2,
      pi_lo = exp(fit$pi_lo), pi_hi = exp(fit$pi_hi), method = fit$method
    )
  }
}
pooled <- do.call(rbind, pooled_rows)
write.csv(pooled, "score_high_pooled.csv", row.names = FALSE)

cat("===== Primary nCRT score-high sensitivity analysis =====\n")
print(pooled[pooled$analysis == "primary_ncrt9", ], row.names = FALSE, digits = 3)

make_forest <- function(dat, title) {
  fit <- pool_random_effects(dat$logOR, dat$se)
  n <- nrow(dat)
  dat$estimate <- exp(dat$logOR)
  dat$lo <- exp(dat$logOR - 1.96 * dat$se)
  dat$hi <- exp(dat$logOR + 1.96 * dat$se)
  dat$y <- n - seq_len(n) + 1.8
  x_min <- min(0.3, min(dat$lo) * 0.9)
  x_max <- max(dat$hi) * 1.3
  diamond <- data.frame(
    x = exp(c(fit$lo, fit$g, fit$hi, fit$g)), y = c(1, 1.28, 1, 0.72)
  )
  ggplot(dat) +
    annotate("rect", xmin = x_min, xmax = x_max, ymin = 0.4, ymax = 1.5,
             fill = "grey95", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi, y = y), orientation = "y",
                  width = 0.18, color = "grey40") +
    geom_point(aes(x = estimate, y = y), shape = 15, size = 2.4, color = "#0072B2") +
    geom_polygon(data = diamond, aes(x = x, y = y), fill = "#D55E00") +
    scale_y_continuous(breaks = c(dat$y, 1), labels = c(dat$cohort, "Pooled (REML-mKH)")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(x_min, x_max)) +
    labs(
      title = sprintf("%s: OR=%.2f (%.2f-%.2f)",
                      title, exp(fit$g), exp(fit$lo), exp(fit$hi)),
      x = "OR (top-tertile score vs rest, log scale)", y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
}

primary <- results[results$treatment_scope == "nCRT", ]
plots <- lapply(features, function(feature) {
  make_forest(primary[primary$feature == feature, ], feature)
})
score_high_grob <- arrangeGrob(grobs = plots, ncol = 3)
png("fig_S_score_high.png", width = 3300, height = 1300, res = 300)
grid::grid.draw(score_high_grob)
dev.off()
ggsave("fig_S_score_high.pdf", score_high_grob, width = 11, height = 4.3,
       device = grDevices::cairo_pdf)

cat("Outputs: immune_hot_xcohort.csv / score_high_pooled.csv / fig_S_score_high.png/.pdf\n")
cat("===== Complete =====\n")
