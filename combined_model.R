# ==============================
# 联合预测模型: IFN + IGF2 + L1CAM (+NK/Cyto)
# 训练: GSE209746 ; 验证: GSE35452, GSE87211
# 所有特征在各队列内 z-score 标准化
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(DESeq2)
  library(pROC)
})

# ---- 基因集 ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]

mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）

# 构建每个队列的特征矩阵: 行=样本, 列=IFN/IGF2/L1CAM/NK/Cyto(全部 z-score)
# symbols_full 保留完整注释(可能含 ///), 基因提取用 token 匹配
build_features <- function(y_raw, symbols_full, grp) {
  sym_first <- sapply(strsplit(symbols_full, " */// *|, *"), `[`, 1)
  y_sym <- symbolize(y_raw, sym_first)
  sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
  sc_mcp <- gsva(ssgseaParam(y_sym, mcp_sets), verbose = FALSE)
  z <- function(x) as.numeric(scale(x))
  get_gene <- function(g) {
    hit <- which(sapply(strsplit(symbols_full, " */// *|, *"),
                        function(x) g %in% x))
    if (length(hit) == 0) stop(g, " 无探针")
    v <- apply(y_raw[hit, , drop = FALSE], 1, var)
    as.numeric(y_raw[hit[which.max(v)], ])
  }
  feat <- data.frame(
    IFN  = z(colMeans(sc_ifn)[colnames(y_raw)]),
    IGF2 = z(get_gene("IGF2")),
    L1CAM = z(get_gene("L1CAM")),
    NK   = z(sc_mcp["NK cells", colnames(y_raw)]),
    Cyto = z(sc_mcp["Cytotoxic lymphocytes", colnames(y_raw)]),
    response = grp[colnames(y_raw)],
    row.names = colnames(y_raw)
  )
  feat
}

# ---- 队列1: GSE209746 (训练) ----
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym1 <- ens2sym[sub("\\..*$", "", rownames(assay(vsd)))]
y1 <- symbolize(assay(vsd), sym1)
grp1 <- ifelse(labels$response[match(colnames(y1), labels$patient)] == "CR", "good", "poor")
names(grp1) <- colnames(y1)
f1 <- build_features(assay(vsd), sym1, grp1)

# ---- 队列2: GSE35452 ----
eset2 <- readRDS("GSE35452_eset.rds")
fd2 <- fData(eset2)
clin2 <- pData(eset2)
rt <- as.character(clin2$characteristics_ch1.1)
grp2 <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
names(grp2) <- colnames(exprs(eset2))
f2 <- build_features(exprs(eset2), as.character(fd2[["Gene symbol"]]), grp2)

# ---- 队列3: GSE87211 ----
eset3 <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
y3raw <- exprs(eset3)
clin3 <- pData(eset3)
tissue_col <- grep("^tissue", colnames(clin3), value = TRUE)[1]
ypT_col <- grep("depth of invasion after", colnames(clin3), value = TRUE)[1]
ypN_col <- grep("lymph node metastasis after", colnames(clin3), value = TRUE)[1]
is_tum <- clin3[[tissue_col]] == "rectal tumor"
ypT <- suppressWarnings(as.numeric(as.character(clin3[[ypT_col]])))
ypN <- suppressWarnings(as.numeric(as.character(clin3[[ypN_col]])))
grp3_all <- ifelse(ypT == 0 & ypN == 0, "good",
                   ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0),
                          "poor", NA))
ok3 <- is_tum & !is.na(grp3_all)
grp3 <- grp3_all[ok3]
stopifnot(sum(grp3 == "good") == 35L, sum(grp3 == "poor") == 167L)
names(grp3) <- colnames(y3raw)[ok3]
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym3 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y3raw)]
f3 <- build_features(y3raw[, ok3], sym3, grp3)

# ==============================
# 训练模型(GSE209746), 结局: poor=1
# ==============================
f1$y <- ifelse(f1$response == "poor", 1, 0)
m_ifn   <- glm(y ~ IFN, data = f1, family = binomial)
m_full  <- glm(y ~ IFN + IGF2 + L1CAM, data = f1, family = binomial)
m_full2 <- glm(y ~ IFN + IGF2 + L1CAM + NK + Cyto, data = f1, family = binomial)

cat("===== 模型系数(GSE209746 训练) =====\n")
print(summary(m_full)$coefficients)

# 落盘三个模型的完整系数表，供论文与独立复核引用（仅记录，不改变任何拟合结果）
coef_to_df <- function(model, model_name) {
  cc <- summary(model)$coefficients
  data.frame(model = model_name, term = rownames(cc),
             estimate = cc[, 1], std_error = cc[, 2],
             statistic = cc[, 3], p_value = cc[, 4],
             row.names = NULL, check.names = FALSE)
}
coef_tbl <- rbind(
  coef_to_df(m_ifn, "IFN"),
  coef_to_df(m_full, "IFN+IGF2+L1CAM"),
  coef_to_df(m_full2, "IFN+IGF2+L1CAM+NK+Cyto"))
write.csv(coef_tbl, "combined_model_coefficients.csv", row.names = FALSE)

# ==============================
# 各队列 AUC 评估
# ==============================
AUC_LOG <- list()
eval_auc <- function(feat, model, label) {
  lp <- predict(model, newdata = feat)
  r <- roc(ifelse(feat$response == "poor", 1, 0), as.numeric(lp),
           levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- ci.auc(r)
  cat(sprintf("%-38s %s AUC = %.3f (%.3f-%.3f)\n",
              label, "", as.numeric(auc(r)), ci[1], ci[3]))
  AUC_LOG[[length(AUC_LOG) + 1L]] <<- data.frame(
    model = CURRENT_MODEL, cohort = trimws(label),
    auc = as.numeric(auc(r)), ci_lo = ci[1], ci_hi = ci[3])
  r
}

cat("\n===== AUC: IFN 单独 vs 联合模型 =====\n")
for (nm in c("IFN单独", "IFN+IGF2+L1CAM", "+NK+Cyto")) {
  model <- switch(nm,
                  "IFN单独" = m_ifn,
                  "IFN+IGF2+L1CAM" = m_full,
                  "+NK+Cyto" = m_full2)
  CURRENT_MODEL <- nm
  cat("---", nm, "---\n")
  eval_auc(f1, model, "  GSE209746(训练)")
  eval_auc(f2, model, "  GSE35452(验证)")
  eval_auc(f3, model, "  GSE87211(验证)")
}
auc_tbl <- do.call(rbind, AUC_LOG)
write.csv(auc_tbl, "combined_model_auc.csv", row.names = FALSE)

# ==============================
# ROC 图: 联合模型在三个队列
# ==============================
png("combined_model_ROC.png", width = 1100, height = 900, res = 130)
plot(roc(ifelse(f1$response=="poor",1,0), predict(m_full, f1), quiet=TRUE),
     col = "black", lwd = 2.5, legacy.axes = TRUE,
     main = "Combined model: IFN + IGF2 + L1CAM")
plot(roc(ifelse(f2$response=="poor",1,0), predict(m_full, f2), quiet=TRUE),
     col = "gray50", lwd = 2.5, add = TRUE)
plot(roc(ifelse(f3$response=="poor",1,0), predict(m_full, f3), quiet=TRUE),
     col = "tomato", lwd = 2.5, add = TRUE)
abline(0, 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = c(sprintf("GSE209746 train  AUC=%.2f",
                          as.numeric(auc(roc(ifelse(f1$response=="poor",1,0),
                                             predict(m_full,f1), quiet=TRUE)))),
                  sprintf("GSE35452 valid  AUC=%.2f",
                          as.numeric(auc(roc(ifelse(f2$response=="poor",1,0),
                                             predict(m_full,f2), quiet=TRUE)))),
                  sprintf("GSE87211 valid  AUC=%.2f",
                          as.numeric(auc(roc(ifelse(f3$response=="poor",1,0),
                                             predict(m_full,f3), quiet=TRUE))))),
       col = c("black", "gray50", "tomato"), lwd = 2.5, cex = 0.9)
dev.off()

# Keep model-specific data separate from the canonical cohort feature tables.
# In particular, `y` is a training-only outcome and must not leak into the
# feature-table schema consumed by validation and meta-analysis.
write.csv(f1, "combined_model_features_GSE209746.csv")
write.csv(f2, "combined_model_features_GSE35452.csv")
write.csv(f3, "combined_model_features_GSE87211.csv")

f1_canonical <- f1
f1_canonical$y <- NULL
write.csv(f1_canonical, "features_GSE209746.csv")
write.csv(f2, "features_GSE35452.csv")
write.csv(f3, "features_GSE87211.csv")
cat("\n===== 完成 =====\n")
