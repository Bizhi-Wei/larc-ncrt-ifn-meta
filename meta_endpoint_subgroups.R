# ==============================
# Response-definition sensitivity (primary nine nCRT cohorts only)
#
# Endpoint classes (locked 2026-09-12 with mentor review):
#   complete_response_endpoints
#     GSE209746 sustained 2-year CR vs iCR
#     GSE87211   pCR = ypT0N0
#   broader_response_classification_endpoints
#     remaining seven nCRT cohorts (TRG-based or sensitivity/resistance)
#
# Output:
#   meta_endpoint_subgroups.csv
#   meta_endpoint_q_between.csv
#   meta_endpoint_effects_classified.csv
#
# Interpretation target is direction and magnitude under each endpoint class,
# not formal equivalence between subgroups. The complete-response subgroup has
# only k = 2: pooled estimate is exploratory; prediction interval is omitted
# because t-based PI uses df = k - 2 = 0.
# ==============================

source("project_config.R", encoding = "UTF-8")
source("meta_utils.R", encoding = "UTF-8")

endpoint_class <- c(
  GSE209746 = "complete_response",
  GSE87211  = "complete_response",
  GSE35452  = "broader_response",
  GSE150082 = "broader_response",
  GSE45404  = "broader_response",
  GSE119409 = "broader_response",
  GSE133057 = "broader_response",
  GSE53781  = "broader_response",
  GSE94104  = "broader_response"
)
class_label <- c(
  complete_response = "complete-response endpoints (sustained CR / pCR)",
  broader_response = "broader response-classification endpoints (TRG or sensitivity/resistance)"
)

features <- c("IFN", "NK", "Cyto")
primary <- LARC_CONFIG$analysis_sets$primary_ncrt9
stopifnot(setequal(primary, names(endpoint_class)))
stopifnot(sum(endpoint_class == "complete_response") == 2L)
stopifnot(sum(endpoint_class == "broader_response") == 7L)

expected_n <- setNames(LARC_CONFIG$cohorts$n, LARC_CONFIG$cohorts$cohort)
n_complete <- sum(expected_n[names(endpoint_class)[endpoint_class == "complete_response"]])
n_broader <- sum(expected_n[names(endpoint_class)[endpoint_class == "broader_response"]])
stopifnot(n_complete == 295L, n_broader == 282L)

# Prefer the already-computed cohort effects when present; otherwise recompute.
if (file.exists("meta_all_cohorts.csv")) {
  effects <- read.csv("meta_all_cohorts.csv", stringsAsFactors = FALSE)
} else {
  effects <- NULL
}

rows <- list()
for (cohort in primary) {
  if (!is.null(effects)) {
    sub <- effects[effects$cohort == cohort & effects$feature %in% features, ]
    stopifnot(nrow(sub) == length(features))
    for (i in seq_len(nrow(sub))) {
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = cohort,
        endpoint_class = endpoint_class[[cohort]],
        feature = sub$feature[i],
        n = sub$n[i],
        n_good = sub$n_good[i],
        n_poor = sub$n_poor[i],
        g = sub$g[i],
        se = sub$se[i],
        stringsAsFactors = FALSE
      )
    }
  } else {
    dat <- read.csv(paste0("features_", cohort, ".csv"), stringsAsFactors = FALSE)
    for (feature in features) {
      eff <- hedges_g(dat[[feature]], dat$response)
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = cohort,
        endpoint_class = endpoint_class[[cohort]],
        feature = feature,
        n = nrow(dat),
        n_good = unname(eff["n_good"]),
        n_poor = unname(eff["n_poor"]),
        g = unname(eff["g"]),
        se = unname(eff["se"]),
        stringsAsFactors = FALSE
      )
    }
  }
}
classified <- do.call(rbind, rows)
classified$lo <- classified$g - 1.96 * classified$se
classified$hi <- classified$g + 1.96 * classified$se
classified <- classified[order(classified$endpoint_class, classified$feature, classified$cohort), ]
write.csv(classified, "meta_endpoint_effects_classified.csv", row.names = FALSE)

pooled_rows <- list()
for (feature in features) {
  for (cls in c("complete_response", "broader_response")) {
    dat <- classified[classified$feature == feature & classified$endpoint_class == cls, ]
    want_pi <- cls == "broader_response"
    fit <- pool_random_effects(dat$g, dat$se, prediction = want_pi)
    pooled_rows[[length(pooled_rows) + 1L]] <- data.frame(
      feature = feature,
      endpoint_class = cls,
      endpoint_class_label = class_label[[cls]],
      k = fit$k,
      n_total = sum(dat$n),
      estimate = fit$g,
      ci_lo = fit$lo,
      ci_hi = fit$hi,
      p = fit$p,
      I2 = fit$I2,
      tau2 = fit$tau2,
      Q = fit$Q,
      Qp = fit$Qp,
      pi_lo = fit$pi_lo,
      pi_hi = fit$pi_hi,
      pi_reported = want_pi,
      direction_positive = sum(dat$g > 0),
      direction_total = nrow(dat),
      method = fit$method,
      interpretation = ifelse(
        cls == "complete_response",
        "exploratory: k=2; interpret for direction/magnitude; no prediction interval",
        "primary-method pooling within broader endpoint class"
      ),
      stringsAsFactors = FALSE
    )
  }
}
pooled <- do.call(rbind, pooled_rows)
write.csv(pooled, "meta_endpoint_subgroups.csv", row.names = FALSE)

# Subgroup-difference test via Cochran Q decomposition (fixed-effect weights).
# A non-significant Q_between is NOT evidence of no difference; the comparison
# is underpowered (especially k=2 complete-response arm).
qb_rows <- list()
for (feature in features) {
  dat <- classified[classified$feature == feature, ]
  w <- 1 / dat$se^2
  mu_fixed <- sum(w * dat$g) / sum(w)
  Q_total <- sum(w * (dat$g - mu_fixed)^2)

  Q_within <- 0
  df_within <- 0
  for (cls in unique(dat$endpoint_class)) {
    s <- dat[dat$endpoint_class == cls, ]
    ws <- 1 / s$se^2
    mus <- sum(ws * s$g) / sum(ws)
    Q_within <- Q_within + sum(ws * (s$g - mus)^2)
    df_within <- df_within + (nrow(s) - 1L)
  }
  Q_between <- Q_total - Q_within
  df_between <- length(unique(dat$endpoint_class)) - 1L
  p_between <- stats::pchisq(max(Q_between, 0), df = df_between, lower.tail = FALSE)

  qb_rows[[length(qb_rows) + 1L]] <- data.frame(
    feature = feature,
    Q_total = Q_total,
    Q_within = Q_within,
    Q_between = max(Q_between, 0),
    df_between = df_between,
    p_between = p_between,
    note = "exploratory; underpowered; do not interpret NS as equivalence"
  )
}
q_between <- do.call(rbind, qb_rows)
write.csv(q_between, "meta_endpoint_q_between.csv", row.names = FALSE)

cat("\n===== Endpoint-class cohort effects (IFN shown fully) =====\n")
print(classified[classified$feature == "IFN", ], row.names = FALSE, digits = 3)
cat("\n===== Pooled by endpoint class =====\n")
print(pooled[, c(
  "feature", "endpoint_class", "k", "n_total", "estimate", "ci_lo", "ci_hi",
  "p", "I2", "tau2", "pi_lo", "pi_hi", "pi_reported", "direction_positive"
)], row.names = FALSE, digits = 3)
cat("\n===== Subgroup-difference Q_between =====\n")
print(q_between, row.names = FALSE, digits = 3)
cat("\nComplete-response subgroup: k=2, no prediction interval by design.\n")
