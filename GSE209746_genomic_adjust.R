# ==============================
# GSE209746 基因组学混杂校正:
#  IFN~response 校正 WES 纯度 / TMB / fraction_cna / genome_doubled
#  (S1C 为论文官方 WES 逐患者数据; 纯度为 WES 推断的真纯度, 强于间质代理)
# ==============================

feat <- read.csv("features_GSE209746.csv", stringsAsFactors = FALSE)
wes  <- read.csv("paper_S1C_wes.csv", stringsAsFactors = FALSE)
names(feat)[1] <- "Patient_ID"
# S1C 的 Patient_ID 形如 "P-00JM195", features 用 "JM195" —— 提取 JM 号对齐
wes$Patient_ID <- sub(".*(JM[0-9]+).*", "\\1", wes$Patient_ID)

d <- merge(feat, wes, by = "Patient_ID")
d <- d[!is.na(d$TMB) & !is.na(d$purity), ]
d$resp <- as.integer(d$response == "good")
d$ifn_z <- as.numeric(scale(d$IFN))
d$tmb_z <- as.numeric(scale(as.numeric(d$TMB)))
d$pur_z <- as.numeric(scale(as.numeric(d$purity)))
d$cna_z <- as.numeric(scale(as.numeric(d$fraction_cna)))
cat("合并后样本:", nrow(d), " CR:", sum(d$resp), "\n")

# TMB / purity 本身与疗效的关系
cat("\n-- 基因组特征与疗效 --\n")
cov_models <- list()
for (v in c("tmb_z", "pur_z", "cna_z")) {
  g <- glm(reformulate(v, "resp"), data = d, family = binomial)
  cov_models[[v]] <- g
  cat(sprintf("%-7s 单变量 P = %.3f\n", v, summary(g)$coefficients[2, 4]))
}
gd <- glm(resp ~ factor(genome_doubled), data = d, family = binomial)
cat(sprintf("genome_doubled 单变量 P = %.3f\n", summary(gd)$coefficients[2, 4]))

# 线性模型: IFN ~ response + 基因组协变量
m0 <- lm(ifn_z ~ resp, data = d)
m1 <- lm(ifn_z ~ resp + pur_z + tmb_z + cna_z, data = d)
cat("\n-- IFN ~ response --\n")
cat(sprintf("未校正:        beta=%.3f P=%.3f\n",
            coef(m0)["resp"], summary(m0)$coefficients["resp", 4]))
cat(sprintf("校正纯度+TMB+CNA: beta=%.3f P=%.3f\n",
            coef(m1)["resp"], summary(m1)$coefficients["resp", 4]))

# 逻辑模型: response ~ IFN + 基因组协变量
g0 <- glm(resp ~ ifn_z, data = d, family = binomial)
g1 <- glm(resp ~ ifn_z + pur_z + tmb_z + cna_z, data = d, family = binomial)
cat("\n-- response ~ IFN --\n")
cat(sprintf("未校正:          OR=%.2f P=%.3f\n",
            exp(coef(g0)["ifn_z"]), summary(g0)$coefficients["ifn_z", 4]))
cat(sprintf("校正纯度+TMB+CNA: OR=%.2f P=%.3f\n",
            exp(coef(g1)["ifn_z"]), summary(g1)$coefficients["ifn_z", 4]))

# MSI 三分层外的终极模型: response ~ IFN + MSI + 纯度 + TMB
d$msi <- read.csv("paper_94_labels.csv")$msi_status[match(d$Patient_ID,
         read.csv("paper_94_labels.csv")$patient)]
g2 <- glm(resp ~ ifn_z + msi + pur_z + tmb_z, data = d[!is.na(d$msi), ], family = binomial)
cat(sprintf("\n-- response ~ IFN + MSI + 纯度 + TMB (n=%d) --\n", nrow(d[!is.na(d$msi), ])))
print(summary(g2)$coefficients)

out <- data.frame(
  model = c("lm_IFN~resp", "lm_IFN~resp+purity+TMB+CNA",
            "glm_resp~IFN", "glm_resp~IFN+purity+TMB+CNA",
            "glm_resp~IFN+MSI+purity+TMB"),
  term = "IFN/response",
  estimate = c(coef(m0)["resp"], coef(m1)["resp"],
               exp(coef(g0)["ifn_z"]), exp(coef(g1)["ifn_z"]),
               exp(coef(g2)["ifn_z"])),
  p = c(summary(m0)$coefficients["resp", 4], summary(m1)$coefficients["resp", 4],
        summary(g0)$coefficients["ifn_z", 4], summary(g1)$coefficients["ifn_z", 4],
        summary(g2)$coefficients["ifn_z", 4]))
write.csv(out, "GSE209746_genomic_adjust.csv", row.names = FALSE)

model_row <- function(fit, model_name, focal_term, effect_measure, adjustment,
                      outcome = "Good response") {
  term_name <- if (focal_term == "Genome doubled") {
    rownames(coef(summary(fit)))[2]
  } else if (focal_term == "Response (good vs poor)") {
    "resp"
  } else {
    "ifn_z"
  }
  if (focal_term == "TMB (per SD)") term_name <- "tmb_z"
  if (focal_term == "WES purity (per SD)") term_name <- "pur_z"
  if (focal_term == "Fraction CNA (per SD)") term_name <- "cna_z"
  est <- unname(coef(fit)[term_name])
  ci <- unname(confint.default(fit)[term_name, ])
  if (effect_measure == "Odds ratio") {
    est <- exp(est)
    ci <- exp(ci)
  }
  mf <- model.frame(fit)
  event_n <- if (outcome == "Good response") sum(model.response(mf) == 1) else sum(d$resp)
  data.frame(
    section = "WES genomic adjustment", cohort = "GSE209746",
    population = "WES-available pretreatment tumors", outcome = outcome,
    model = model_name, focal_term = focal_term, effect_measure = effect_measure,
    estimate = est, ci_lo = ci[1], ci_hi = ci[2],
    p_value = coef(summary(fit))[term_name, 4], n = nobs(fit), events = event_n,
    adjustment = adjustment,
    note = ifelse(effect_measure == "Odds ratio",
                  "OR above 1 indicates greater odds of good response.",
                  "Positive beta indicates higher IFN score in responders."),
    source = "GSE209746_genomic_adjust_models.csv",
    stringsAsFactors = FALSE
  )
}

gd_term <- rownames(coef(summary(gd)))[2]
genomic_detail <- rbind(
  model_row(cov_models[["tmb_z"]], "Univariable logistic regression", "TMB (per SD)",
            "Odds ratio", "None"),
  model_row(cov_models[["pur_z"]], "Univariable logistic regression", "WES purity (per SD)",
            "Odds ratio", "None"),
  model_row(cov_models[["cna_z"]], "Univariable logistic regression", "Fraction CNA (per SD)",
            "Odds ratio", "None"),
  model_row(gd, "Univariable logistic regression", "Genome doubled",
            "Odds ratio", "None"),
  model_row(m0, "Linear regression", "Response (good vs poor)", "Beta", "None",
            outcome = "IFN score (z)"),
  model_row(m1, "Linear regression", "Response (good vs poor)", "Beta",
            "WES purity, TMB, and fraction CNA", outcome = "IFN score (z)"),
  model_row(g0, "Logistic regression", "IFN score (per SD)", "Odds ratio", "None"),
  model_row(g1, "Logistic regression", "IFN score (per SD)", "Odds ratio",
            "WES purity, TMB, and fraction CNA"),
  model_row(g2, "Logistic regression", "IFN score (per SD)", "Odds ratio",
            "MSI status, WES purity, and TMB")
)
write.csv(genomic_detail, "GSE209746_genomic_adjust_models.csv", row.names = FALSE)
cat("\n输出: GSE209746_genomic_adjust.csv / GSE209746_genomic_adjust_models.csv\n===== 完成 =====\n")
