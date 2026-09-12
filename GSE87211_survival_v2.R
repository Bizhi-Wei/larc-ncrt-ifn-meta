# ==============================
# GSE87211 生存分析 v2: 多变量Cox / 非线性检验(样条) / pCR分层
# DFS 主结局, OS 补充; IFN/NK/Cyto z分数
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(survival)
  library(splines)
})

eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset)
clin <- pData(eset)

tissue_col <- grep("^tissue", colnames(clin), value = TRUE)[1]
ypT_col <- grep("depth of invasion after", colnames(clin), value = TRUE)[1]
ypN_col <- grep("lymph node metastasis after", colnames(clin), value = TRUE)[1]
is_tum <- clin[[tissue_col]] == "rectal tumor"
# 临床随访列含空字符串(""), as.numeric 会发 "NAs introduced by coercion";
# 空值随后由 !is.na(...) & time>0 剔除, 故显式抑制该已知良性提示。
num_clin <- function(x) suppressWarnings(as.numeric(as.character(x)))
dfs_time  <- num_clin(clin[["disease free time (month):ch1"]])
dfs_event <- num_clin(clin[["cancer recurrance after surgery:ch1"]])
os_time   <- num_clin(clin[["survival time (month):ch1"]])
os_event  <- num_clin(clin[["death due to tumor:ch1"]])
ypT <- suppressWarnings(as.numeric(as.character(clin[[ypT_col]])))
ypN <- suppressWarnings(as.numeric(as.character(clin[[ypN_col]])))
age <- num_clin(gsub("^age: ", "", as.character(clin[["age:ch1"]])))
gender <- as.character(clin[["gender:ch1"]])

ok <- is_tum & !is.na(dfs_time) & !is.na(dfs_event) & dfs_time > 0
cat("DFS 可用肿瘤样本:", sum(ok), " 事件:", sum(dfs_event[ok]), "\n")

# ---- 表达特征 ----
y_tum <- y[, ok]
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym_full <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y_tum)]
sym_first <- sapply(strsplit(sym_full, " */// *|, *"), `[`, 1)
source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）
y_sym <- symbolize(y_tum, sym_first)

gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
sc_mcp <- gsva(ssgseaParam(y_sym, mcp_sets), verbose = FALSE)

srv <- data.frame(
  time = dfs_time[ok], event = dfs_event[ok],
  ifn = as.numeric(colMeans(sc_ifn)),
  nk  = as.numeric(sc_mcp["NK cells", ]),
  cyto = as.numeric(sc_mcp["Cytotoxic lymphocytes", ]),
  ypT = ypT[ok], ypN = ypN[ok],
  age = age[ok], gender = gender[ok],
  os_time = os_time[ok], os_event = os_event[ok]
)
srv$pCR <- ifelse(srv$ypT == 0 & srv$ypN == 0, "pCR",
           ifelse((!is.na(srv$ypT) & srv$ypT != 0) |
                    (!is.na(srv$ypN) & srv$ypN != 0), "non-pCR", NA))
srv$ypN_bin <- factor(ifelse(srv$ypN >= 1, "Npos", "N0"))
srv$ypN_bin <- relevel(srv$ypN_bin, ref = "N0")
srv$ifn_z <- as.numeric(scale(srv$ifn))
srv$nk_z <- as.numeric(scale(srv$nk))
srv$cyto_z <- as.numeric(scale(srv$cyto))

phr <- function(fit, label) {
  s <- summary(fit)
  for (v in rownames(s$coefficients)) {
    cat(sprintf("  %-28s %-12s HR = %.3f (%.3f-%.3f)  P = %.4f\n",
                label, v, s$conf.int[v, 1], s$conf.int[v, 3], s$conf.int[v, 4],
                s$coefficients[v, 5]))
  }
}

cox_row <- function(fit, endpoint, model_name, predictor, term_name, adjustment,
                    population = "Pretreatment nCRT tumors with follow-up") {
  s <- summary(fit)
  data.frame(
    cohort = "GSE87211", endpoint = endpoint, population = population,
    model = model_name, predictor = predictor, effect_measure = "Hazard ratio per SD",
    estimate = unname(s$conf.int[term_name, 1]),
    ci_lo = unname(s$conf.int[term_name, 3]), ci_hi = unname(s$conf.int[term_name, 4]),
    p_value = unname(s$coefficients[term_name, 5]), n = fit$n, events = fit$nevent,
    adjustment = adjustment, note = "HR below 1 indicates lower event hazard.",
    source = "GSE87211_survival_models.csv", stringsAsFactors = FALSE
  )
}

# ==============================
# 1. 多变量 Cox (DFS): 校正 yp 分期
# ==============================
cat("\n===== 多变量 Cox (DFS) =====\n")
mv <- srv[!is.na(srv$ypT) & !is.na(srv$ypN), ]
cat("多变量模型样本:", nrow(mv), " 事件:", sum(mv$event), "\n")
cat("-- 模型A: IFN + ypT(有序) + ypN(N+ vs N0) --\n")
fitA <- coxph(Surv(time, event) ~ ifn_z + ypT + ypN_bin, data = mv)
phr(fitA, "A")
cat("-- 模型B: A + 年龄 + 性别 --\n")
fitB <- coxph(Surv(time, event) ~ ifn_z + ypT + ypN_bin + age + gender, data = mv)
phr(fitB, "B")
cat("-- 模型C: NK/Cyto 同样校正 --\n")
fitC1 <- coxph(Surv(time, event) ~ nk_z + ypT + ypN_bin, data = mv)
phr(fitC1, "C-nk")
fitC2 <- coxph(Surv(time, event) ~ cyto_z + ypT + ypN_bin, data = mv)
phr(fitC2, "C-cyto")

# ==============================
# 2. 非线性检验: 自然样条 vs 线性
# ==============================
cat("\n===== IFN 非线性检验(ns df=3 vs 线性, DFS) =====\n")
fit_lin <- coxph(Surv(time, event) ~ ifn_z, data = srv)
fit_ns  <- coxph(Surv(time, event) ~ ns(ifn_z, df = 3), data = srv)
lr <- anova(fit_lin, fit_ns)
lr_p <- lr[2, grep("Chi|Pr", colnames(lr), value = TRUE)]
lr_p <- as.numeric(lr_p[!is.na(lr_p)][length(lr_p[!is.na(lr_p)])])
cat(sprintf("似然比检验 P = %.4f (>0.05 支持线性假设)\n", lr_p))

# 样条曲线图: 相对中位数的 HR
png("GSE87211_IFN_DFS_spline.png", width = 1000, height = 850, res = 130)
grid <- seq(min(srv$ifn_z), max(srv$ifn_z), length.out = 100)
nsb <- ns(srv$ifn_z, df = 3)
grid_ns <- predict(nsb, grid)
coefs <- coef(fit_ns)
lp <- as.numeric(grid_ns %*% coefs)
ref <- as.numeric(predict(nsb, 0) %*% coefs)   # 以 z=0(均值)为参照
hr <- exp(lp - ref)
plot(grid, hr, type = "l", lwd = 2.5, col = "steelblue",
     xlab = "IFN score (z)", ylab = "Hazard ratio (vs mean)",
     main = sprintf("GSE87211 DFS: IFN dose-response (nonlinearity P = %.3f)",
                    lr_p))
abline(h = 1, lty = 2, col = "gray50")
rug(srv$ifn_z, col = "gray70")
dev.off()

# ==============================
# 3. pCR 分层: 各层内 IFN 预后价值
# ==============================
cat("\n===== 按 pCR 分层(DFS) =====\n")
fit_nonpcr <- NULL
nonpcr_n <- nonpcr_events <- NA_integer_
pcr_n <- pcr_events <- NA_integer_
for (st in c("pCR", "non-pCR")) {
  sub <- srv[srv$pCR == st, ]
  if (st == "pCR") {
    pcr_n <- nrow(sub); pcr_events <- sum(sub$event)
  } else {
    nonpcr_n <- nrow(sub); nonpcr_events <- sum(sub$event)
  }
  cat(sprintf("-- %s: n=%d, 事件=%d --\n", st, nrow(sub), sum(sub$event)))
  if (sum(sub$event) >= 5) {
    fit <- coxph(Surv(time, event) ~ ifn_z, data = sub)
    if (st == "non-pCR") fit_nonpcr <- fit
    phr(fit, st)
    med <- ifelse(sub$ifn > median(sub$ifn), "high", "low")
    sd2 <- survdiff(Surv(time, event) ~ med, data = sub)
    cat(sprintf("   KM中位切分 log-rank P = %.4f\n", 1 - pchisq(sd2$chisq, 1)))
  } else {
    cat("   事件数<5, 仅描述: 高IFN组事件",
        sum(sub$event[sub$ifn > median(sub$ifn)]), "/",
        sum(sub$ifn > median(sub$ifn)), " 低IFN组事件",
        sum(sub$event[sub$ifn <= median(sub$ifn)]), "/",
        sum(sub$ifn <= median(sub$ifn)), "\n")
  }
}

# ==============================
# 4. OS 多变量补充
# ==============================
cat("\n===== OS 多变量 =====\n")
ok_os <- !is.na(srv$os_time) & !is.na(srv$os_event) & srv$os_time > 0 &
         !is.na(srv$ypT) & !is.na(srv$ypN)
mo <- srv[ok_os, ]
cat("OS 样本:", nrow(mo), " 事件:", sum(mo$os_event), "\n")
fit_os <- coxph(Surv(os_time, os_event) ~ ifn_z + ypT + ypN_bin, data = mo)
phr(fit_os, "OS-A")

# ---- 汇总保存 ----
summ <- data.frame(
  model = c("DFS_univariate_IFN", "DFS_multivar_A_IFN", "DFS_multivar_B_IFN",
            "DFS_nonlinearity_LR", "DFS_OS_multivar_IFN"),
  stat = c(
    sprintf("HR=%.3f P=%.4f", summary(fit_lin)$conf.int[1], summary(fit_lin)$coefficients[5]),
    sprintf("HR=%.3f P=%.4f", summary(fitA)$conf.int["ifn_z", 1],
            summary(fitA)$coefficients["ifn_z", 5]),
    sprintf("HR=%.3f P=%.4f", summary(fitB)$conf.int["ifn_z", 1],
            summary(fitB)$coefficients["ifn_z", 5]),
    sprintf("P=%.4f", lr_p),
    sprintf("HR=%.3f P=%.4f", summary(fit_os)$conf.int["ifn_z", 1],
            summary(fit_os)$coefficients["ifn_z", 5]))
)
write.csv(summ, "GSE87211_survival_v2_summary.csv", row.names = FALSE)

survival_detail <- rbind(
  cox_row(fit_lin, "DFS", "Univariable Cox", "IFN score", "ifn_z", "None"),
  cox_row(fitA, "DFS", "Multivariable Cox A", "IFN score", "ifn_z",
          "ypT (ordinal) and ypN (positive vs N0)"),
  cox_row(fitB, "DFS", "Multivariable Cox B", "IFN score", "ifn_z",
          "ypT, ypN, age, and sex"),
  cox_row(fitC1, "DFS", "Multivariable Cox", "NK score", "nk_z",
          "ypT (ordinal) and ypN (positive vs N0)"),
  cox_row(fitC2, "DFS", "Multivariable Cox", "Cytotoxic score", "cyto_z",
          "ypT (ordinal) and ypN (positive vs N0)"),
  data.frame(
    cohort = "GSE87211", endpoint = "DFS",
    population = "Pretreatment nCRT tumors with follow-up",
    model = "Likelihood-ratio test", predictor = "IFN natural spline vs linear term",
    effect_measure = "P value only", estimate = NA_real_, ci_lo = NA_real_,
    ci_hi = NA_real_, p_value = lr_p, n = nrow(srv), events = sum(srv$event),
    adjustment = "None", note = "Tests non-linearity using a 3-df natural spline.",
    source = "GSE87211_survival_models.csv", stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE87211", endpoint = "DFS", population = "pCR tumors",
    model = "Cox model not estimated", predictor = "IFN score",
    effect_measure = "Hazard ratio per SD", estimate = NA_real_, ci_lo = NA_real_,
    ci_hi = NA_real_, p_value = NA_real_, n = pcr_n, events = pcr_events,
    adjustment = "None", note = "Not estimated because no DFS events occurred in pCR tumors.",
    source = "GSE87211_survival_models.csv", stringsAsFactors = FALSE
  ),
  cox_row(fit_nonpcr, "DFS", "pCR-stratified univariable Cox", "IFN score", "ifn_z",
          "None", population = "Non-pCR tumors"),
  cox_row(fit_os, "OS", "Multivariable Cox", "IFN score", "ifn_z",
          "ypT (ordinal) and ypN (positive vs N0)")
)
write.csv(survival_detail, "GSE87211_survival_models.csv", row.names = FALSE)

cat("\n输出: GSE87211_survival_v2_summary.csv / GSE87211_survival_models.csv\n===== 完成 =====\n")
