# ==============================
# GSE87211: pCR vs non-pCR limma 分析
# + 论文15基因方向交叉验证
# pCR 定义: ypT0 且 ypN0 (CR类似) ; non-pCR (iCR类似)
# logFC > 0 = pCR 高表达 => 论文 iCR 高表达基因在此应 logFC < 0
# ==============================

suppressMessages({
  library(GEOquery)
  library(limma)
})

# ---- 1. 读取数据 ----
eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
expr <- exprs(eset)
clin <- pData(eset)

cat("表达矩阵维度:", dim(expr), "\n")
cat("range:", paste(round(range(expr), 2), collapse = " ~ "), "\n")

# ---- 2. 肿瘤样本 + pCR 分组 ----
tum <- clin$tissue_ch1 == "rectal tumor" |
       (("tissue:ch1" %in% colnames(clin)) & FALSE)
# 列名兼容处理
tissue_col <- grep("^tissue", colnames(clin), value = TRUE)[1]
ypT_col <- grep("depth of invasion after", colnames(clin), value = TRUE)[1]
ypN_col <- grep("lymph node metastasis after", colnames(clin), value = TRUE)[1]

is_tum <- clin[[tissue_col]] == "rectal tumor"
ypT <- suppressWarnings(as.numeric(as.character(clin[[ypT_col]])))
ypN <- suppressWarnings(as.numeric(as.character(clin[[ypN_col]])))
# 三态判定：任一明确非0即可判定 non-pCR；只有 ypT0N0 才是 pCR；其余缺失。
group_all <- ifelse(ypT == 0 & ypN == 0, "pCR",
                    ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0),
                           "non-pCR", NA))
ok <- is_tum & !is.na(group_all)
group <- group_all[ok]
stopifnot(sum(group == "pCR") == 35L, sum(group == "non-pCR") == 167L)
cat("\n===== 分组 =====\n")
print(table(group))

expr_tum <- expr[, ok]
group <- factor(group, levels = c("non-pCR", "pCR"))

# ---- 3. 如需要则 log2 ----
if (max(expr_tum) > 100) {
  cat("最大值>100, 做 log2 转换\n")
  expr_tum <- log2(expr_tum + 1)
}

# ---- 4. limma ----
design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr_tum, design))
res <- topTable(fit, coef = "grouppCR", number = Inf, sort.by = "none")

cat("\n===== 差异分析概况(pCR vs non-pCR) =====\n")
cat("P<0.05:", sum(res$P.Value < 0.05),
    " adj.P<0.05:", sum(res$adj.P.Val < 0.05), "\n")

# ---- 5. 探针注释 ----
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym_vec <- setNames(pmap$GENE_SYMBOL, pmap$ID)
res$symbol <- sym_vec[rownames(res)]

# ---- 6. 15基因方向交叉验证 ----
genes15 <- data.frame(
  symbol = c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300"),
  paper_dir = c(rep("iCR_high", 8), rep("CR_high", 7)),
  stringsAsFactors = FALSE)

cat("\n===== 15基因在 GSE87211 中的方向 =====\n")
cat(sprintf("%-9s %-16s %8s %9s   %-9s %s\n",
            "symbol", "best_probe", "logFC", "P.Value", "本队列方向", "一致?"))
n_same <- 0; n_opp <- 0; n_na <- 0
out <- data.frame()
for (i in seq_len(nrow(genes15))) {
  g <- genes15$symbol[i]
  syms <- strsplit(res$symbol, " */// *| *, *")
  hit <- rownames(res)[sapply(syms, function(x) g %in% x)]
  if (length(hit) == 0) {
    cat(sprintf("%-9s %-16s %8s %9s   无探针\n", g, "NA", "NA", "NA"))
    n_na <- n_na + 1
    out <- rbind(out, data.frame(symbol = g, best_probe = NA, logFC = NA,
                                 P.Value = NA, dir_here = NA,
                                 paper_dir = genes15$paper_dir[i],
                                 match = "no_probe"))
    next
  }
  sub <- res[hit, ]
  best <- sub[which.max(abs(sub$t)), ]
  d <- ifelse(best$logFC > 0, "CR_high", "iCR_high")
  same <- ifelse(d == genes15$paper_dir[i], "same", "opposite")
  if (same == "same") n_same <- n_same + 1 else n_opp <- n_opp + 1
  cat(sprintf("%-9s %-16s %+8.3f %9.4f   %-9s %s\n",
              g, rownames(best), best$logFC, best$P.Value, d, same))
  out <- rbind(out, data.frame(symbol = g, best_probe = rownames(best),
                               logFC = round(best$logFC, 3),
                               P.Value = signif(best$P.Value, 3),
                               dir_here = d,
                               paper_dir = genes15$paper_dir[i],
                               match = same))
}
cat("\n方向一致:", n_same, " 相反:", n_opp, " 无探针:", n_na, "\n")

# ---- 7. IGF2 / L1CAM 全部探针 ----
cat("\n===== IGF2 / L1CAM 全部探针 =====\n")
for (g in c("IGF2", "L1CAM")) {
  syms <- strsplit(res$symbol, " */// *| *, *")
  hit <- rownames(res)[sapply(syms, function(x) g %in% x)]
  cat("---", g, "---\n")
  if (length(hit) > 0) print(round(res[hit, c("logFC","AveExpr","t","P.Value")], 3))
  else cat("无探针\n")
}

write.csv(res, "GSE87211_limma_pCR_vs_nonpCR.csv")
write.csv(out, "GSE87211_15genes_crosscheck.csv", row.names = FALSE)
saveRDS(eset, "GSE87211_eset.rds")
cat("\n===== 完成 =====\n")
