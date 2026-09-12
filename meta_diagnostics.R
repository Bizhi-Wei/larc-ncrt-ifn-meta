# Diagnostics for the primary nine-cohort nCRT meta-analysis.
# Produces leave-one-cohort-out estimates, Egger tests, and correctly centred funnel plots.

suppressMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})
source("project_config.R", encoding = "UTF-8")
source("meta_utils.R")

effects <- read.csv("meta_all_cohorts.csv", stringsAsFactors = FALSE)
effects <- effects[effects$cohort %in% LARC_CONFIG$analysis_sets$primary_ncrt9, ]
features <- c("IFN", "NK", "Cyto")
stopifnot(setequal(unique(effects$cohort), LARC_CONFIG$analysis_sets$primary_ncrt9))

loo_rows <- list()
for (feature in features) {
  dat <- effects[effects$feature == feature, ]
  full <- pool_random_effects(dat$g, dat$se)
  for (i in seq_len(nrow(dat))) {
    fit <- pool_random_effects(dat$g[-i], dat$se[-i])
    loo_rows[[length(loo_rows) + 1L]] <- data.frame(
      feature = feature, omitted = dat$cohort[i], k = fit$k,
      g = fit$g, lo = fit$lo, hi = fit$hi, p = fit$p,
      I2 = fit$I2, tau2 = fit$tau2, method = fit$method
    )
  }
  loo_rows[[length(loo_rows) + 1L]] <- data.frame(
    feature = feature, omitted = "(none, full model)", k = full$k,
    g = full$g, lo = full$lo, hi = full$hi, p = full$p,
    I2 = full$I2, tau2 = full$tau2, method = full$method
  )
}
loo <- do.call(rbind, loo_rows)
write.csv(loo, "meta_leave_one_out.csv", row.names = FALSE)

egger_rows <- list()
for (feature in features) {
  dat <- effects[effects$feature == feature, ]
  standardized <- dat$g / dat$se
  precision <- 1 / dat$se
  fit <- stats::lm(standardized ~ precision)
  coefficient <- summary(fit)$coefficients
  egger_rows[[feature]] <- data.frame(
    feature = feature, k = nrow(dat),
    egger_intercept = coefficient[1, 1],
    egger_t = coefficient[1, 3],
    egger_p = coefficient[1, 4],
    note = sprintf("k=%d; low power; interpret cautiously", nrow(dat))
  )
}
egger <- do.call(rbind, egger_rows)
rownames(egger) <- NULL
write.csv(egger, "meta_egger.csv", row.names = FALSE)

full_rows <- loo[loo$omitted == "(none, full model)", c("feature", "g")]
colnames(full_rows)[2] <- "full_g"
loo_plot_data <- merge(
  loo[loo$omitted != "(none, full model)", ], full_rows, by = "feature"
)
p_loo <- ggplot(loo_plot_data, aes(x = g, y = reorder(omitted, g))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  geom_vline(aes(xintercept = full_g), color = "#D55E00", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15, color = "grey40") +
  geom_point(shape = 15, color = "#0072B2", size = 2) +
  facet_wrap(~feature, scales = "free_y", nrow = 1) +
  labs(x = "Pooled Hedges g after omitting one cohort", y = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))

make_funnel <- function(dat, title) {
  fit <- pool_random_effects(dat$g, dat$se)
  se_sequence <- seq(0.01, max(dat$se) * 1.15, length.out = 100)
  funnel <- data.frame(
    se = se_sequence,
    lo = fit$g - 1.96 * se_sequence,
    hi = fit$g + 1.96 * se_sequence
  )
  ggplot(dat, aes(x = g, y = se)) +
    geom_line(data = funnel, aes(x = lo, y = se), linetype = 2, color = "grey60") +
    geom_line(data = funnel, aes(x = hi, y = se), linetype = 2, color = "grey60") +
    geom_vline(xintercept = fit$g, color = "#D55E00", linewidth = 0.5) +
    geom_point(size = 2.3, color = "#0072B2") +
    geom_text(aes(label = cohort), vjust = -0.75, size = 2.2, check_overlap = TRUE) +
    scale_y_reverse() +
    labs(title = title, x = "Hedges g", y = "Standard error") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"))
}

funnel_plots <- lapply(features, function(feature) {
  make_funnel(effects[effects$feature == feature, ], feature)
})
funnel_row <- arrangeGrob(grobs = funnel_plots, ncol = 3)

png("fig_S_funnel.png", width = 3000, height = 1100, res = 300)
grid.draw(funnel_row)
dev.off()

png("fig_S_meta_diagnostics.png", width = 3000, height = 2200, res = 300)
grid.arrange(p_loo, funnel_row, ncol = 1, heights = c(1.1, 1))
dev.off()
diagnostic_grob <- arrangeGrob(p_loo, funnel_row, ncol = 1, heights = c(1.1, 1))
ggsave("fig_S_meta_diagnostics.pdf", diagnostic_grob, width = 10, height = 7.3,
       device = grDevices::cairo_pdf)

cat("Outputs: meta_leave_one_out.csv / meta_egger.csv / ",
    "fig_S_funnel.png / fig_S_meta_diagnostics.png/.pdf\n", sep = "")
cat("===== Complete =====\n")
