# ==============================
# TCGA-READ: IFN 评分预后价值 (非新辅助队列, 旁证)
# 数据: Xena GDC star_counts (原始count) + Liu2018 curated survival + clinical matrix
# 注意: ssGSEA 为样本内秩次方法, count 的单调变换不影响结果
# ==============================

suppressMessages({
  library(GSVA)
  library(survival)
})

# ---- 表达: star counts ----
cnt <- read.delim("TCGA-READ.star_counts.tsv.gz", row.names = 1, check.names = FALSE)
cat("基因数:", nrow(cnt), " 样本数:", ncol(cnt), "\n")
# 仅原发肿瘤 (barcode 第14-15位 = 01)
sam_type <- substr(colnames(cnt), 14, 15)
cnt <- cnt[, sam_type == "01"]
cat("原发肿瘤样本:", ncol(cnt), "\n")

# Ensembl → symbol (复用论文 S8A 注释)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym <- ens2sym[sub("\\..*$", "", rownames(cnt))]
cat("Ensembl→symbol 注释覆盖率:", round(mean(!is.na(sym)) * 100, 1), "%\n")

y <- log2(as.matrix(cnt) + 1)
source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）
y_sym <- symbolize(y, sym)
cat("symbol 化后基因数:", nrow(y_sym), "\n")

# ---- 评分 ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
sc_mcp <- gsva(ssgseaParam(y_sym, mcp_sets), verbose = FALSE)

srv_df <- data.frame(
  sample = colnames(y_sym),
  ifn = as.numeric(colMeans(sc_ifn)),
  nk  = as.numeric(sc_mcp["NK cells", colnames(y_sym)]),
  cyto = as.numeric(sc_mcp["Cytotoxic lymphocytes", colnames(y_sym)])
)

# ---- 生存与临床 ----
surv <- read.delim("READ_survival.txt", stringsAsFactors = FALSE, check.names = FALSE)
surv <- surv[!duplicated(surv$`_PATIENT`), ]   # 每患者一行(弃 -11 正常重复行)
clin <- read.delim("READ_clinicalMatrix", row.names = 1, check.names = FALSE)
srv_df$patient <- substr(srv_df$sample, 1, 12)
srv_df <- merge(srv_df, surv[, c("_PATIENT", "OS", "OS.time", "PFI", "PFI.time")],
                by.x = "patient", by.y = "_PATIENT")
srv_df$stage <- clin[srv_df$sample, "pathologic_stage"]
srv_df$neoadj <- clin[srv_df$sample, "history_of_neoadjuvant_treatment"]
srv_df$rad <- clin[srv_df$sample, "additional_radiation_therapy"]
if (all(is.na(srv_df$stage))) {   # 样本名不匹配时退回患者级匹配
  srv_df$stage <- clin[match(srv_df$patient, substr(rownames(clin), 1, 12)), "pathologic_stage"]
  srv_df$neoadj <- clin[match(srv_df$patient, substr(rownames(clin), 1, 12)), "history_of_neoadjuvant_treatment"]
  cat("!! 样本名未匹配, 已退回患者级匹配\n")
}

cat("\n合并后样本:", nrow(srv_df), "\n")
cat("新辅助治疗史:", paste(names(table(srv_df$neoadj, useNA = "ifany")),
                          table(srv_df$neoadj, useNA = "ifany"), collapse = " / "), "\n")
cat("分期分布:", paste(names(table(srv_df$stage, useNA = "ifany")),
                      table(srv_df$stage, useNA = "ifany"), collapse = " / "), "\n")

# 排除新辅助治疗史(保持治疗初治原发队列), NA 视为未治疗并做敏感性
srv_main <- srv_df[srv_df$neoadj != "Yes" | is.na(srv_df$neoadj), ]
# 分期 → 有序数值
st <- tolower(trimws(srv_main$stage))
srv_main$stage_n <- NA
srv_main$stage_n[st %in% c("stage i")] <- 1
srv_main$stage_n[st %in% c("stage ii", "stage iia", "stage iib", "stage iic")] <- 2
srv_main$stage_n[st %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic")] <- 3
srv_main$stage_n[st %in% c("stage iv", "stage iva", "stage ivb")] <- 4
cat("主分析样本(排除新辅助):", nrow(srv_main),
    " OS事件:", sum(srv_main$OS == 1, na.rm = TRUE),
    " PFI事件:", sum(srv_main$PFI == 1, na.rm = TRUE), "\n")

srv_main$ifn_z <- as.numeric(scale(srv_main$ifn))
srv_main$nk_z <- as.numeric(scale(srv_main$nk))
srv_main$cyto_z <- as.numeric(scale(srv_main$cyto))

phr <- function(fit, label) {
  s <- summary(fit)
  for (v in rownames(s$coefficients)) {
    cat(sprintf("  %-18s %-8s HR = %.3f (%.3f-%.3f)  P = %.4f\n",
                label, v, s$conf.int[v, 1], s$conf.int[v, 3], s$conf.int[v, 4],
                s$coefficients[v, 5]))
  }
}

cox_row <- function(fit, endpoint, model_name, predictor, term_name, adjustment) {
  s <- summary(fit)
  data.frame(
    cohort = "TCGA-READ", endpoint = endpoint,
    population = "Primary tumors without documented neoadjuvant treatment",
    model = model_name, predictor = predictor, effect_measure = "Hazard ratio per SD",
    estimate = unname(s$conf.int[term_name, 1]),
    ci_lo = unname(s$conf.int[term_name, 3]), ci_hi = unname(s$conf.int[term_name, 4]),
    p_value = unname(s$coefficients[term_name, 5]), n = fit$n, events = fit$nevent,
    adjustment = adjustment,
    note = "Treatment-naive prognostic context; not an nCRT validation cohort.",
    source = "TCGA_READ_survival_models.csv", stringsAsFactors = FALSE
  )
}

# ---- OS / PFI: 单变量 + 校正分期 ----
model_rows <- list()
row_i <- 0L
for (endpoint in c("OS", "PFI")) {
  cat(sprintf("\n===== %s =====\n", endpoint))
  tt <- (if (endpoint == "OS") srv_main$OS.time else srv_main$PFI.time) / 30.44
  ee <- if (endpoint == "OS") srv_main$OS else srv_main$PFI
  ok <- !is.na(tt) & !is.na(ee) & tt > 0
  d <- srv_main[ok, ]; d$tt <- tt[ok]; d$ee <- ee[ok]
  cat("可用:", nrow(d), " 事件:", sum(d$ee), "\n")
  for (v in c("ifn_z", "nk_z", "cyto_z")) {
    fml <- reformulate(v, response = "Surv(tt, ee)")
    fit <- coxph(fml, data = d)
    s <- summary(fit)
    cat(sprintf("  单变量 %-7s HR = %.3f (%.3f-%.3f)  P = %.4f\n",
                v, s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5]))
    row_i <- row_i + 1L
    predictor_label <- c(ifn_z = "IFN score", nk_z = "NK score", cyto_z = "Cytotoxic score")[[v]]
    model_rows[[row_i]] <- cox_row(fit, endpoint, "Univariable Cox", predictor_label,
                                    v, "None")
  }
  d2 <- d[!is.na(d$stage_n), ]
  cat("  -- 校正分期(n=", nrow(d2), ") --\n")
  fit_mv <- coxph(Surv(tt, ee) ~ ifn_z + stage_n, data = d2)
  phr(fit_mv, "多变量")
  row_i <- row_i + 1L
  model_rows[[row_i]] <- cox_row(fit_mv, endpoint, "Multivariable Cox", "IFN score",
                                  "ifn_z", "AJCC stage (ordinal)")
}

# ---- KM: IFN 中位切分 (OS) ----
ok <- !is.na(srv_main$OS.time) & !is.na(srv_main$OS) & srv_main$OS.time > 0
d <- srv_main[ok, ]
d$os_m <- d$OS.time / 30.44
d$ifn_grp <- ifelse(d$ifn > median(d$ifn), "IFN-high", "IFN-low")
sdiff <- survdiff(Surv(os_m, OS) ~ ifn_grp, data = d)
p_km <- 1 - pchisq(sdiff$chisq, 1)
cat(sprintf("\nOS KM 中位切分: log-rank P = %.4f\n", p_km))
print(table(d$ifn_grp, d$OS))

png("TCGA_READ_IFN_OS_KM.png", width = 1000, height = 850, res = 130)
km <- survfit(Surv(os_m, OS) ~ ifn_grp, data = d)
plot(km, col = c("tomato", "steelblue"), lwd = 2.5,
     xlab = "Months", ylab = "Overall survival",
     main = sprintf("TCGA-READ OS by IFN score (log-rank P = %.3f)", p_km),
     mark.time = TRUE)
legend("bottomleft", legend = c("IFN-high", "IFN-low"),
       col = c("tomato", "steelblue"), lwd = 2.5)
dev.off()

write.csv(srv_main, "TCGA_READ_ifn_survival.csv", row.names = FALSE)
write.csv(do.call(rbind, model_rows), "TCGA_READ_survival_models.csv", row.names = FALSE)
cat("\n输出: TCGA_READ_ifn_survival.csv / TCGA_READ_survival_models.csv\n===== 完成 =====\n")
