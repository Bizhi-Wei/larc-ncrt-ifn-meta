# 宽松过滤复跑: 总count >= 10(接近论文做法,依靠DESeq2独立过滤)
suppressMessages(library(DESeq2))

exp_raw <- as.matrix(read.csv(gzfile("GSE209746_geneExp_n_114_ensembl.csv.gz"),
                              row.names = 1, check.names = FALSE))
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)

count_mat <- exp_raw[, labels$patient]
count_filt <- count_mat[rowSums(count_mat) >= 10, ]
cat("过滤后基因数:", nrow(count_filt), "\n")

coldata <- data.frame(
  response = factor(labels$response, levels = c("CR", "iCR")),
  row.names = labels$patient
)
dds <- DESeqDataSetFromMatrix(count_filt, coldata, design = ~ response)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("response", "iCR", "CR"), alpha = 0.05))
res$ensembl_gene_id <- rownames(res)
res$symbol <- ens2sym[sub("\\..*$", "", res$ensembl_gene_id)]

sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
cat("显著基因数(|log2FC|>1 & padj<0.05):", nrow(sig), "\n\n")

genes15 <- c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300")
paper15 <- s8a[s8a$hgnc_symbol %in% genes15, ]
n_rep <- sum(paper15$ensembl_gene_id %in% sig$ensembl_gene_id)
cat("复现个数:", n_rep, "/15\n\n")

for (g in genes15) {
  o <- res[!is.na(res$symbol) & res$symbol == g, ]
  p <- s8a[!is.na(s8a$hgnc_symbol) & s8a$hgnc_symbol == g, ]
  tag <- ifelse(!is.na(o$padj[1]) && o$padj[1] < 0.05 &&
                abs(o$log2FoldChange[1]) > 1, "[复现]", "[---]")
  cat(sprintf("%-9s 我们: log2FC=%+.2f padj=%-12.2e | 论文: log2FC=%+.2f padj=%-12.2e %s\n",
              g, o$log2FoldChange[1], o$padj[1],
              as.numeric(p$log2FoldChange[1]), as.numeric(p$padj[1]), tag))
}

# 整体方向一致性(全部基因)
both <- merge(res, s8a, by = "ensembl_gene_id", suffixes = c("_ours", "_paper"))
both <- both[!is.na(both$padj_ours) & both$padj_paper != "", ]
both$padj_paper <- as.numeric(both$padj_paper)
both$lfc_paper <- as.numeric(both$log2FoldChange_paper)
paper_sig <- both[is.finite(both$lfc_paper) & is.finite(both$padj_paper) &
                  abs(both$lfc_paper) > 1 & both$padj_paper < 0.05, ]
same <- sum(sign(paper_sig$log2FoldChange_ours) == sign(paper_sig$lfc_paper), na.rm = TRUE)
cat("\n按论文阈值匹配到", nrow(paper_sig), "个基因, 我们结果方向一致:",
    same, "/", nrow(paper_sig), "\n")

write.csv(res, "GSE209746_DESeq2_iCR_vs_CR_loose_filter.csv", row.names = FALSE)
saveRDS(dds, "GSE209746_dds_loose.rds")
cat("\n===== 完成 =====\n")
