# =====================================================================
# Figure S2: primary nine-cohort meta-analysis diagnostics
# A: leave-one-cohort-out estimates
# B: centred funnel plots
# C: Baujat heterogeneity contribution and influence
#
# Core conclusion: no single cohort overturns the primary IFN association;
# influential cohorts remain identifiable and are reported transparently.
# Outputs: fig_S2_diagnostics.png/.pdf/.svg
# =====================================================================

suppressMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})
source("project_config.R", encoding = "UTF-8")
source("meta_utils.R")

features <- c("IFN", "Cyto", "NK")
feature_labels <- c(IFN = "IFN score", Cyto = "Cytotoxic lymphocytes", NK = "NK cells")
cohort_order <- LARC_CONFIG$analysis_sets$primary_ncrt9

effects <- read.csv("meta_all_cohorts.csv", stringsAsFactors = FALSE)
effects <- effects[effects$cohort %in% cohort_order & effects$feature %in% features, ]
loo <- read.csv("meta_leave_one_out.csv", stringsAsFactors = FALSE)
baujat <- read.csv("meta_baujat.csv", stringsAsFactors = FALSE)
stopifnot(setequal(unique(effects$cohort), cohort_order))

full_rows <- loo[loo$omitted == "(none, full model)", c("feature", "g")]
colnames(full_rows)[2] <- "full_g"
loo_df <- merge(loo[loo$omitted != "(none, full model)", ], full_rows, by = "feature")
loo_df <- loo_df[loo_df$feature %in% features, ]
loo_df$feature <- factor(loo_df$feature, levels = features, labels = feature_labels[features])
loo_df$omitted <- factor(loo_df$omitted, levels = rev(cohort_order))

p_loo <- ggplot(loo_df, aes(x = g, y = omitted)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70", linewidth = 0.35) +
  geom_vline(aes(xintercept = full_g), color = "#D55E00", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15,
                color = "grey40", linewidth = 0.35) +
  geom_point(shape = 15, color = "#0072B2", size = 1.7) +
  facet_wrap(~feature, nrow = 1) +
  labs(x = "Pooled Hedges g after omitting one cohort", y = NULL) +
  theme_classic(base_size = 7.5) +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.text = element_text(color = "black"))

make_funnel <- function(feature) {
  dat <- effects[effects$feature == feature, ]
  fit <- pool_random_effects(dat$g, dat$se)
  se_sequence <- seq(0.01, max(dat$se) * 1.15, length.out = 100)
  funnel <- data.frame(
    se = se_sequence,
    lo = fit$g - 1.96 * se_sequence,
    hi = fit$g + 1.96 * se_sequence
  )
  ggplot(dat, aes(x = g, y = se)) +
    geom_line(data = funnel, aes(x = lo, y = se), linetype = 2,
              color = "grey60", linewidth = 0.35) +
    geom_line(data = funnel, aes(x = hi, y = se), linetype = 2,
              color = "grey60", linewidth = 0.35) +
    geom_vline(xintercept = fit$g, color = "#D55E00", linewidth = 0.45) +
    geom_point(size = 1.8, color = "#0072B2") +
    geom_text(aes(label = cohort), vjust = -0.75, size = 2.0, check_overlap = TRUE) +
    scale_y_reverse(expand = expansion(mult = c(0.04, 0.16))) +
    labs(title = feature_labels[[feature]], x = "Hedges g", y = "Standard error") +
    theme_classic(base_size = 7.5) +
    theme(plot.title = element_text(face = "bold", size = 8),
          axis.text = element_text(color = "black"))
}
funnel_row <- arrangeGrob(grobs = lapply(features, make_funnel), ncol = 3)

make_baujat <- function(feature) {
  dat <- baujat[baujat$feature == feature, ]
  label_idx <- union(which.max(dat$Q_contrib), which.max(dat$influence))
  dat$label <- ""
  dat$label[label_idx] <- dat$cohort[label_idx]
  eff <- effects[effects$feature == feature, ]
  fit <- pool_random_effects(eff$g, eff$se)
  ggplot(dat, aes(x = Q_contrib, y = influence)) +
    geom_point(size = 1.8, color = "#0072B2") +
    geom_text(data = dat[dat$label != "", ], aes(label = label),
              vjust = -0.8, size = 2.1, check_overlap = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(title = sprintf("%s\nQ=%.2f, tau2=%.3f, I2=%.0f%%",
                         feature_labels[[feature]], fit$Q, fit$tau2, fit$I2),
         x = "Contribution to heterogeneity (Q)",
         y = "Influence on pooled estimate") +
    theme_classic(base_size = 7.5) +
    theme(plot.title = element_text(face = "bold", size = 7.4, lineheight = 0.95),
          axis.text = element_text(color = "black"))
}
baujat_row <- arrangeGrob(grobs = lapply(features, make_baujat), ncol = 3)

row_title <- function(label) textGrob(label, x = 0, hjust = 0,
                                      gp = gpar(fontface = "bold", fontsize = 9))
full_grob <- arrangeGrob(
  arrangeGrob(p_loo, top = row_title("A  Leave-one-cohort-out sensitivity")),
  arrangeGrob(funnel_row, top = row_title("B  Centred funnel plots")),
  arrangeGrob(baujat_row, top = row_title("C  Baujat influence diagnostics")),
  ncol = 1, heights = c(1.1, 1, 1)
)

width_in <- 183 / 25.4
height_in <- 240 / 25.4
png("fig_S2_diagnostics.png", width = width_in, height = height_in,
    units = "in", res = 300)
grid.draw(full_grob)
dev.off()
grDevices::cairo_pdf("fig_S2_diagnostics.pdf", width = width_in, height = height_in,
                     family = "Arial")
grid.draw(full_grob)
dev.off()
svglite::svglite("fig_S2_diagnostics.svg", width = width_in, height = height_in)
grid.draw(full_grob)
dev.off()

cat("Outputs: fig_S2_diagnostics.png/.pdf/.svg\n")
cat("===== Complete =====\n")
