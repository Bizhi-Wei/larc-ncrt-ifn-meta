# Consolidate prediction intervals produced by the primary REML-mKH analyses.
# Run meta_all_cohorts.R and immune_hot_xcohort.R first.

score <- read.csv("meta_pooled_results.csv", stringsAsFactors = FALSE)
score_out <- data.frame(
  analysis = "score_meta", feature = score$feature, subset = score$analysis,
  k = score$k, n_total = score$n_total, est = score$estimate,
  ci_lo = score$ci_lo, ci_hi = score$ci_hi, p = score$p,
  I2 = score$I2, tau2 = score$tau2,
  pi_lo = score$pi_lo, pi_hi = score$pi_hi, method = score$method
)

score_high <- read.csv("score_high_pooled.csv", stringsAsFactors = FALSE)
score_high_out <- data.frame(
  analysis = "score_high_OR", feature = score_high$feature, subset = score_high$analysis,
  k = score_high$k, n_total = score_high$n_total, est = score_high$OR,
  ci_lo = score_high$ci_lo, ci_hi = score_high$ci_hi, p = score_high$p,
  I2 = score_high$I2, tau2 = score_high$tau2_log,
  pi_lo = score_high$pi_lo, pi_hi = score_high$pi_hi, method = score_high$method
)

out <- rbind(score_out, score_high_out)
write.csv(out, "meta_prediction_intervals.csv", row.names = FALSE)
print(out, row.names = FALSE, digits = 3)
cat("Output: meta_prediction_intervals.csv\n===== Complete =====\n")
