# ==============================
# GSE35452 差异表达分析
# Responder vs Non-responder
# ==============================

library(GEOquery)
library(limma)

eset <- readRDS("GSE35452_eset.rds")
expr_matrix <- exprs(eset)
clinical <- pData(eset)

# ---- 分组(和上次一样) ----
response_text <- as.character(clinical$characteristics_ch1.1)
group <- ifelse(
  grepl("Non-responder", response_text, ignore.case = TRUE),
  "Non-responder",
  "Responder"
)
group <- factor(group, levels = c("Non-responder", "Responder"))
cat("分组人数:\n")
print(table(group))

# ---- ① range ----
cat("\n===== ① range(expr_matrix) =====\n")
print(range(expr_matrix))

# ---- ② 箱线图 ----
png("GSE35452_boxplot.png", width = 1200, height = 800, res = 120)
boxplot(expr_matrix,
        las = 2, outline = FALSE,
        col = ifelse(group == "Responder", "tomato", "steelblue"),
        main = "GSE35452 Expression Boxplot",
        ylab = "Expression value")
legend("topright",
       legend = c("Non-responder", "Responder"),
       fill = c("steelblue", "tomato"))
dev.off()

# ---- ③ PCA 图 ----
pca <- prcomp(t(expr_matrix), scale. = TRUE)
var_pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

png("GSE35452_PCA.png", width = 900, height = 800, res = 120)
plot(pca$x[, 1], pca$x[, 2],
     col = ifelse(group == "Responder", "tomato", "steelblue"),
     pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", var_pct[1], "%)"),
     ylab = paste0("PC2 (", var_pct[2], "%)"),
     main = "GSE35452 PCA")
legend("topright",
       legend = c("Non-responder", "Responder"),
       col = c("steelblue", "tomato"),
       pch = 19)
dev.off()

# ---- ④ limma 差异分析 ----
design <- model.matrix(~ group)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "groupResponder",
                number = Inf, sort.by = "none")

n_005   <- sum(res$adj.P.Val < 0.05)
n_005_05 <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5)
n_005_1 <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)

cat("\n===== ④ 差异基因数量 =====\n")
cat("adj.P.Val < 0.05                    :", n_005, "\n")
cat("adj.P.Val < 0.05 & |logFC| > 0.5    :", n_005_05, "\n")
cat("adj.P.Val < 0.05 & |logFC| > 1      :", n_005_1, "\n")

cat("\nP.Value < 0.05 未校正的数量(参考):", sum(res$P.Value < 0.05), "\n")
cat("最小 adj.P.Val:", min(res$adj.P.Val), "\n")

write.csv(res, "GSE35452_limma_results.csv")

cat("\n===== 完成 =====\n")
