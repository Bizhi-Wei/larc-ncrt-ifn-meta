# ==============================
# GSE209746: CR vs iCR DESeq2 分析
# 复现 Nature Medicine 2022 (PMID 35970919)
# 标签采用论文 Table S1A 定义(Use_for_Response_2year=1)
# ==============================

suppressMessages({
  library(DESeq2)
})

# ---- 1. 读取 count 矩阵 ----
exp_raw <- as.matrix(
  read.csv(
    gzfile("GSE209746_geneExp_n_114_ensembl.csv.gz"),
    row.names = 1, check.names = FALSE
  )
)

# ---- 2. 读论文标签(93例) + 论文注释(S8A ensembl->symbol) ----
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
cat("论文标签例数:", nrow(labels), "\n")
print(table(labels$response))

s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)

# ---- 3. 建立 count matrix(93例,列序与标签一致) ----
common <- intersect(labels$patient, colnames(exp_raw))
cat("能匹配表达矩阵的患者:", length(common), "\n")
labels <- labels[match(common, labels$patient), ]
count_mat <- exp_raw[, labels$patient]
stopifnot(identical(colnames(count_mat), labels$patient))
cat("count matrix 维度:", dim(count_mat), "\n")

# ---- 4. 过滤低表达基因 ----
min_group <- min(table(labels$response))
keep_gene <- rowSums(count_mat >= 10) >= min_group
count_filt <- count_mat[keep_gene, ]
cat("\n===== 低表达过滤 =====\n")
cat("标准: 至少", min_group, "个样本 count >= 10\n")
cat("过滤前:", nrow(count_mat), " 过滤后:", nrow(count_filt), "\n")

# ---- 5. DESeq2 ----
coldata <- data.frame(
  response = factor(labels$response, levels = c("CR", "iCR")),
  row.names = labels$patient
)

dds <- DESeqDataSetFromMatrix(count_filt, coldata, design = ~ response)
dds <- DESeq(dds)

cat("\n===== size factors =====\n")
print(summary(sizeFactors(dds)))

# ---- 6. VST + PCA ----
vsd <- vst(dds, blind = TRUE)

msi <- labels$msi_status
png("GSE209746_PCA_vst.png", width = 900, height = 800, res = 120)
pd <- plotPCA(vsd, intgroup = "response", returnData = TRUE)
pv <- round(100 * attr(pd, "percentVar"))
plot(pd$PC1, pd$PC2,
     col = ifelse(labels$response == "iCR", "tomato", "steelblue"),
     pch = ifelse(msi == "1", 17, 19),
     cex = 1.4,
     xlab = paste0("PC1 (", pv[1], "%)"),
     ylab = paste0("PC2 (", pv[2], "%)"),
     main = "GSE209746 PCA (VST, 93 samples)")
legend("topright",
       legend = c("iCR", "CR", "MSS(圆)", "MSI(三角)"),
       col = c("tomato", "steelblue", "gray40", "gray40"),
       pch = c(19, 19, 19, 17))
dev.off()

# ---- 7. 差异分析: iCR vs CR(log2FC>0 = iCR高表达,与论文同向) ----
res <- results(dds, contrast = c("response", "iCR", "CR"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
res_df$symbol <- ens2sym[sub("\\..*$", "", res_df$ensembl_gene_id)]
res_df <- res_df[order(res_df$pvalue), ]

cat("\n===== 我们的结果概况 =====\n")
cat("padj < 0.05 基因数:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("padj < 0.05 & |log2FC| > 1:", sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n")

sig_ours <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 &
                   abs(res_df$log2FoldChange) > 1, ]
cat("\n===== 我们的显著基因(|log2FC|>1 & padj<0.05) =====\n")
print(sig_ours[, c("symbol", "ensembl_gene_id", "baseMean",
                   "log2FoldChange", "pvalue", "padj")], row.names = FALSE)

# ---- 8. 对照论文 15 个蛋白编码基因 ----
genes15 <- c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300")

cat("\n===== 论文15基因: 我们 vs 论文(S8A) =====\n")
for (g in genes15) {
  ours <- res_df[res_df$symbol == g, ]
  paper <- s8a[s8a$hgnc_symbol == g, ]
  if (nrow(ours) > 0 && nrow(paper) > 0) {
    cat(sprintf("%-9s 我们: log2FC=%+.2f padj=%.2e | 论文: log2FC=%+.2f padj=%.2e\n",
                g, ours$log2FoldChange[1], ours$padj[1],
                as.numeric(paper$log2FoldChange[1]), as.numeric(paper$padj[1])))
  } else if (nrow(ours) == 0) {
    cat(sprintf("%-9s 我们的矩阵中不存在(可能被过滤或不在GEO的19662基因中)\n", g))
  }
}

# 复现统计
paper15_ens <- s8a$ensembl_gene_id[s8a$hgnc_symbol %in% genes15]
ours_sig_ens <- sig_ours$ensembl_gene_id
n_reproduced <- sum(paper15_ens %in% ours_sig_ens)
cat("\n论文15基因中, 我们也达到 |log2FC|>1 & padj<0.05 的个数:", n_reproduced, "/ 15\n")

# 方向一致性(只看存在于我们结果中的)
both <- merge(res_df, s8a, by = "ensembl_gene_id", suffixes = c("_ours", "_paper"))
both15 <- both[both$symbol %in% genes15, ]
same_dir <- sum(sign(both15$log2FoldChange_ours) == sign(as.numeric(both15$log2FoldChange_paper)))
cat("15基因中存在于我们结果且方向一致:", same_dir, "/", nrow(both15), "\n")

# ---- 9. IGF2 / L1CAM ----
cat("\n===== IGF2 / L1CAM =====\n")
print(res_df[res_df$symbol %in% c("IGF2", "L1CAM"),
             c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj")],
      row.names = FALSE)

# ---- 保存 ----
write.csv(res_df, "GSE209746_DESeq2_iCR_vs_CR.csv", row.names = FALSE)
saveRDS(dds, "GSE209746_dds.rds")
saveRDS(vsd, "GSE209746_vsd.rds")

cat("\n===== 完成 =====\n")
