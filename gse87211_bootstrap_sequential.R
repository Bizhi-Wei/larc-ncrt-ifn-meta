# ==============================
# ⑧ bootstrap 临床模型 + 生存序贯调整 (导师清单⑧)
# A) GSE87211 pCR~cT+cN vs +IFN: 2000 次 bootstrap 优化校正 AUC 差
#    (Efron-Gong optimism bootstrap)
# B) DFS Cox 序贯: ypT,ypN -> +IFN 的 LR 检验与 HR 轨迹
# 输出: gse87211_bootstrap_auc.csv / gse87211_sequential_survival.csv
# ==============================

set.seed(20260912)
suppressMessages({ library(survival); library(Biobase) })
setwd("D:/GEO2R")

# ---- 数据: 复用 clinical_utility 的解析 ----
eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
pd <- pData(eset)
tis <- pd[[grep("^tissue", colnames(pd), value = TRUE)[1]]]
ct_col <- grep("depth of invasion before", colnames(pd), value = TRUE)[1]
cn_col <- grep("lymph node metastasis before", colnames(pd), value = TRUE)[1]
ypT_col <- grep("depth of invasion after", colnames(pd), value = TRUE)[1]
ypN_col <- grep("lymph node metastasis after", colnames(pd), value = TRUE)[1]
ct <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(pd[[ct_col]]))))
cn <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(pd[[cn_col]]))))
ypT <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(pd[[ypT_col]]))))
ypN <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(pd[[ypN_col]]))))
pcr <- ifelse(ypT == 0 & ypN == 0, 1L,
              ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0), 0L, NA))
feat <- read.csv("features_GSE87211.csv", row.names = 1)
ifn <- setNames(feat$IFN[match(colnames(eset), rownames(feat))], colnames(eset))

ok <- tis == "rectal tumor" & !is.na(pcr) & !is.na(ct) & !is.na(cn) & !is.na(ifn)
d <- data.frame(pcr = pcr[ok], ct = ct[ok], cn = cn[ok], ifn = as.numeric(ifn[ok]))
cat("建模样本:", nrow(d), " pCR:", sum(d$pcr), "\n")

auc <- function(y, lp) {
  # Mann-Whitney AUC, 方向: lp 越大越可能是 y=1
  r <- rank(lp); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fit_lp <- function(dat) {
  f1 <- glm(pcr ~ ct + cn, data = dat, family = binomial)
  f2 <- glm(pcr ~ ct + cn + ifn, data = dat, family = binomial)
  list(lp1 = predict(f1, type = "response"), lp2 = predict(f2, type = "response"))
}

# ---- 表观 AUC 与 LR ----
lp <- fit_lp(d)
app1 <- auc(d$pcr, lp$lp1); app2 <- auc(d$pcr, lp$lp2)
f1 <- glm(pcr ~ ct + cn, data = d, family = binomial)
f2 <- glm(pcr ~ ct + cn + ifn, data = d, family = binomial)
lr_p <- anova(f1, f2, test = "LRT")[["Pr(>Chi)"]][2]
cat(sprintf("表观 AUC: cT+cN=%.3f, +IFN=%.3f, delta=%.3f, LR P=%.4g\n",
            app1, app2, app2 - app1, lr_p))

# ---- Optimism bootstrap (2000 次) ----
B <- 2000
opt <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample(nrow(d), replace = TRUE)
  db <- d[idx, ]
  lpb <- tryCatch(fit_lp(db), error = function(e) NULL)
  if (is.null(lpb)) next
  # bootstrap 样本内表现
  inb2 <- auc(db$pcr, lpb$lp2); inb1 <- auc(db$pcr, lpb$lp1)
  # 原样本表现
  lpo <- tryCatch({
    f1b <- glm(pcr ~ ct + cn, data = db, family = binomial)
    f2b <- glm(pcr ~ ct + cn + ifn, data = db, family = binomial)
    list(lp1 = predict(f1b, newdata = d, type = "response"),
         lp2 = predict(f2b, newdata = d, type = "response"))
  }, error = function(e) NULL)
  if (is.null(lpo)) next
  outb2 <- auc(d$pcr, lpo$lp2); outb1 <- auc(d$pcr, lpo$lp1)
  opt[b] <- (inb2 - inb1) - (outb2 - outb1)
}
opt_valid <- opt[is.finite(opt)]
delta_star <- app2 - app1
delta_corr <- delta_star - mean(opt_valid)
boot_se <- sd(quant_opt <- opt_valid)
# percentile CI for corrected delta
delta_boot <- delta_star - opt_valid
ci <- quantile(delta_boot, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("optimism-corrected delta AUC = %.3f (95%% CI %.3f-%.3f), optimism mean=%.4f\n",
            delta_corr, ci[1], ci[2], mean(opt_valid)))
write.csv(data.frame(apparent_base = app1, apparent_full = app2,
                     delta_apparent = delta_star, optimism_mean = mean(opt_valid),
                     delta_corrected = delta_corr, ci_lo = ci[1], ci_hi = ci[2],
                     lr_p = lr_p, B = B, n = nrow(d)),
          "gse87211_bootstrap_auc.csv", row.names = FALSE)

# ---- DFS 生存序贯调整 ----
dfs_t <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(
  pd[["disease free time (month):ch1"]]))))
dfs_e <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(
  pd[["cancer recurrance after surgery:ch1"]]))))
oks <- tis == "rectal tumor" & !is.na(dfs_t) & dfs_t > 0 & !is.na(dfs_e) &
  !is.na(ypT) & !is.na(ypN) & !is.na(ifn)
ds <- data.frame(time = dfs_t[oks], event = dfs_e[oks],
                 ypT = ypT[oks], ypN = ypN[oks], ifn = as.numeric(ifn[oks]))
cat("\n生存样本:", nrow(ds), "事件:", sum(ds$event), "\n")
m1 <- coxph(Surv(time, event) ~ ypT + ypN, data = ds)
m2 <- coxph(Surv(time, event) ~ ypT + ypN + ifn, data = ds)
lrp2 <- anova(m1, m2, test = "LRT")[["Pr(>Chi)"]][2]
# 进一步 + 年龄/性别
age <- suppressWarnings(as.numeric(sub("^[^:]*: *", "", as.character(
  pd[["age:ch1"]]))))
sex <- sub("^[^:]*: *", "", as.character(pd[["gender:ch1"]]))
m3 <- coxph(Surv(time, event) ~ ypT + ypN + ifn + age + sex, data = transform(ds,
  age = age[oks], sex = sex[oks]))
seq <- data.frame(
  model = c("M1: ypT+ypN", "M2: M1+IFN", "M3: M2+age+sex"),
  n = c(nrow(ds), nrow(ds), nrow(ds)),
  events = c(sum(ds$event), sum(ds$event), sum(ds$event)),
  cindex = round(c(summary(m1)$concordance["C"], summary(m2)$concordance["C"],
                   summary(m3)$concordance["C"]), 3),
  ifn_hr = round(c(NA, exp(coef(m2)["ifn"]), exp(coef(m3)["ifn"])), 3),
  ifn_ci = c(NA,
             sprintf("%.2f-%.2f", exp(confint(m2)["ifn", 1]), exp(confint(m2)["ifn", 2])),
             sprintf("%.2f-%.2f", exp(confint(m3)["ifn", 1]), exp(confint(m3)["ifn", 2]))),
  ifn_p = c(NA, round(summary(m2)$coefficients["ifn", "Pr(>|z|)"], 3),
            round(summary(m3)$coefficients["ifn", "Pr(>|z|)"], 3)),
  lr_vs_prev_p = c(NA, lrp2, anova(m2, m3, test = "LRT")[["Pr(>Chi)"]][2]))
write.csv(seq, "gse87211_sequential_survival.csv", row.names = FALSE)
cat("\n===== DFS 序贯调整 =====\n"); print(seq, row.names = FALSE)
cat("\n输出: gse87211_bootstrap_auc.csv / gse87211_sequential_survival.csv\n===== 完成 =====\n")
