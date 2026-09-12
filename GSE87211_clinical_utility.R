# ==============================
# GSE87211 临床增益分析:
#  A) cT/cN ± IFN 预测 pCR: AUC 比较 + 似然比检验 + DeLong
#  B) IFN~pCR 校正 KRAS 突变 (遗传混杂)
#  C) 中位随访时间 (reverse KM)
# ==============================

suppressMessages({
  library(GEOquery); library(GSVA); library(pROC); library(survival)
})

eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset)
clin <- pData(eset)

is_tum <- clin[["tissue:ch1"]] == "rectal tumor"
ypT <- suppressWarnings(as.numeric(as.character(clin[["depth of invasion after rct:ch1"]])))
ypN <- suppressWarnings(as.numeric(as.character(clin[["lymph node metastasis after rct:ch1"]])))
cT  <- suppressWarnings(as.numeric(as.character(clin[["depth of invasion before rct:ch1"]])))
cN  <- suppressWarnings(as.numeric(as.character(clin[["lymph node metastasis before rct:ch1"]])))
cM  <- suppressWarnings(as.numeric(as.character(clin[["metastasis before rct:ch1"]])))
kras_raw <- as.character(clin[["kras mutation:ch1"]])
kras <- ifelse(kras_raw == "WT", 0, ifelse(kras_raw %in% c("NA", ""), NA, 1))
# 随访列含空字符串, 空值随后由 !is.na(...) 剔除; 抑制已知良性 coercion 提示。
num_clin <- function(x) suppressWarnings(as.numeric(as.character(x)))
dfs_time  <- num_clin(clin[["disease free time (month):ch1"]])
dfs_event <- num_clin(clin[["cancer recurrance after surgery:ch1"]])

# ---- IFN 分数(与 survival_v2 完全一致的流程) ----
y_tum <- y[, is_tum]
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
ifn <- colMeans(sc_ifn)

d <- data.frame(sample = colnames(y_tum), ifn = as.numeric(ifn[colnames(y_tum)]),
                ypT = ypT[is_tum], ypN = ypN[is_tum], cT = cT[is_tum],
                cN = cN[is_tum], cM = cM[is_tum], kras = kras[is_tum],
                dfs_time = dfs_time[is_tum], dfs_event = dfs_event[is_tum])
d$pcr <- ifelse(d$ypT == 0 & d$ypN == 0, 1L,
         ifelse((!is.na(d$ypT) & d$ypT != 0) |
                  (!is.na(d$ypN) & d$ypN != 0), 0L, NA_integer_))
d <- d[!is.na(d$pcr), ]
stopifnot(nrow(d) == 202L, sum(d$pcr) == 35L)
d$ifn_z <- as.numeric(scale(d$ifn))
cat("肿瘤样本:", nrow(d), " pCR:", sum(d$pcr), "\n")

# ---- A) 临床增益 ----
da <- d[!is.na(d$cT) & !is.na(d$cN), ]
cat("cTNM 完整样本:", nrow(da), " pCR:", sum(da$pcr), "\n")

m_clin <- glm(pcr ~ cT + cN, data = da, family = binomial)
m_full <- glm(pcr ~ cT + cN + ifn_z, data = da, family = binomial)
m_ifn  <- glm(pcr ~ ifn_z, data = da, family = binomial)

lr <- anova(m_clin, m_full, test = "LRT")
lr_p <- lr[["Pr(>Chi)"]][2]

roc_clin <- roc(da$pcr, predict(m_clin, type = "response"), quiet = TRUE)
roc_full <- roc(da$pcr, predict(m_full, type = "response"), quiet = TRUE)
roc_ifn  <- roc(da$pcr, predict(m_ifn, type = "response"), quiet = TRUE)
delong <- roc.test(roc_clin, roc_full, method = "delong")

cat(sprintf("AUC 临床(cT+cN): %.3f (%.3f-%.3f)\n",
            auc(roc_clin), ci.auc(roc_clin)[1], ci.auc(roc_clin)[3]))
cat(sprintf("AUC IFN 单独:    %.3f (%.3f-%.3f)\n",
            auc(roc_ifn), ci.auc(roc_ifn)[1], ci.auc(roc_ifn)[3]))
cat(sprintf("AUC 临床+IFN:    %.3f (%.3f-%.3f)\n",
            auc(roc_full), ci.auc(roc_full)[1], ci.auc(roc_full)[3]))
cat(sprintf("LR 检验(+IFN): P = %.4g ; DeLong: P = %.4g\n", lr_p, delong$p.value))
cat("全模型系数:\n"); print(summary(m_full)$coefficients)

# ---- B) KRAS 校正 ----
dk <- d[!is.na(d$kras), ]
cat(sprintf("\nKRAS 可用: %d 例 (突变 %d / WT %d), pCR: %d\n",
            nrow(dk), sum(dk$kras == 1), sum(dk$kras == 0), sum(dk$pcr)))
cat("KRAS 突变 vs pCR: ", sprintf("突变 %d/%d, WT %d/%d",
    sum(dk$pcr[dk$kras == 1]), sum(dk$kras == 1),
    sum(dk$pcr[dk$kras == 0]), sum(dk$kras == 0)),
    " Fisher P =", fisher.test(table(dk$kras, dk$pcr))$p.value, "\n")
m_k <- glm(pcr ~ ifn_z + kras, data = dk, family = binomial)
cat("pCR ~ IFN + KRAS:\n"); print(summary(m_k)$coefficients)

# ---- C) 中位随访 (reverse KM, DFS) ----
ds <- d[!is.na(d$dfs_time) & !is.na(d$dfs_event) & d$dfs_time > 0, ]
rk <- survfit(Surv(dfs_time, 1 - dfs_event) ~ 1, data = ds)
med_fu <- summary(rk)$table["median"]
cat(sprintf("\nDFS 队列 n=%d, 事件=%d; 中位随访(reverse KM)= %.1f 月; DFS时间中位= %.1f 月 (range %.1f-%.1f)\n",
            nrow(ds), sum(ds$dfs_event), med_fu, median(ds$dfs_time),
            min(ds$dfs_time), max(ds$dfs_time)))

# ---- 汇总输出 ----
out <- data.frame(item = c("n_tumor", "n_pcr", "n_cTNM", "AUC_clin", "AUC_ifn",
                           "AUC_clin_ifn", "LR_p", "DeLong_p",
                           "KRAS_mut_pcr_rate", "KRAS_wt_pcr_rate", "KRAS_fisher_p",
                           "IFN_in_KRASadj_p", "median_followup_months"),
                  value = c(nrow(d), sum(d$pcr), nrow(da),
                            sprintf("%.3f", auc(roc_clin)), sprintf("%.3f", auc(roc_ifn)),
                            sprintf("%.3f", auc(roc_full)), sprintf("%.4g", lr_p),
                            sprintf("%.4g", delong$p.value),
                            sprintf("%d/%d", sum(dk$pcr[dk$kras == 1]), sum(dk$kras == 1)),
                            sprintf("%d/%d", sum(dk$pcr[dk$kras == 0]), sum(dk$kras == 0)),
                            sprintf("%.4g", fisher.test(table(dk$kras, dk$pcr))$p.value),
                            sprintf("%.4g", summary(m_k)$coefficients["ifn_z", 4]),
                            sprintf("%.1f", med_fu)))
write.csv(out, "GSE87211_clinical_utility.csv", row.names = FALSE)

logit_row <- function(fit, term_name, model_name, adjustment, note) {
  ci <- exp(confint.default(fit)[term_name, ])
  data.frame(
    section = ifelse(grepl("KRAS", adjustment), "KRAS adjustment", "Clinical-stage adjustment"),
    cohort = "GSE87211", population = "Pretreatment tumors with complete model covariates",
    outcome = "Pathological complete response", model = model_name,
    focal_term = ifelse(term_name == "ifn_z", "IFN score (per SD)", "KRAS-mutant vs wild type"),
    effect_measure = "Odds ratio", estimate = unname(exp(coef(fit)[term_name])),
    ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
    p_value = coef(summary(fit))[term_name, 4], n = nobs(fit),
    events = sum(model.response(model.frame(fit)) == 1), adjustment = adjustment,
    note = note, source = "GSE87211_clinical_utility_models.csv",
    stringsAsFactors = FALSE
  )
}
clinical_detail <- rbind(
  logit_row(m_full, "ifn_z", "Logistic regression", "Clinical T and N stage",
            "In-sample incremental association; no optimism correction."),
  logit_row(m_k, "ifn_z", "Logistic regression", "KRAS mutation status",
            "OR above 1 indicates greater odds of pCR per SD higher IFN score."),
  logit_row(m_k, "kras", "Logistic regression", "IFN score and KRAS mutation status",
            "KRAS-mutant tumors are compared with wild-type tumors.")
)
write.csv(clinical_detail, "GSE87211_clinical_utility_models.csv", row.names = FALSE)
cat("\n输出: GSE87211_clinical_utility.csv / GSE87211_clinical_utility_models.csv\n===== 完成 =====\n")
