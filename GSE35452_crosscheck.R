# ==============================
# 交叉验证: 论文15基因在 GSE35452 中的表达方向
# GSE35452: Responder(≈CR) vs Non-responder(≈iCR)
# 方向换算: GSE209746 中 iCR 高表达(log2FC>0)
#   => GSE35452 中应 Non-responder 高表达 => logFC(Res vs NonRes) < 0
# ==============================

suppressMessages(library(Biobase))

eset <- readRDS("GSE35452_eset.rds")
fd <- fData(eset)
limma_res <- read.csv("GSE35452_limma_results.csv", row.names = 1,
                      check.names = FALSE)

# 探针 -> 基因符号
fd$probe <- rownames(fd)
cat("注释用的符号列名:\n")
print(colnames(fd)[1:5])

# 论文15基因(GSE209746 log2FC>0 = iCR高表达)
genes15 <- data.frame(
  symbol = c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300"),
  paper_direction = c(rep("iCR高表达", 8), rep("CR高表达", 7)),
  stringsAsFactors = FALSE
)

sym_col <- "Gene symbol"
results <- list()

for (i in seq_len(nrow(genes15))) {
  g <- genes15$symbol[i]
  # Gene symbol 列可能有 "A /// B" 形式, 精确匹配第一个符号
  syms <- strsplit(as.character(fd[[sym_col]]), " /// ")
  hit <- fd$probe[sapply(syms, function(x) g %in% x)]
  if (length(hit) == 0) {
    results[[g]] <- data.frame(
      symbol = g, n_probe = 0, probes = NA, best_probe = NA,
      logFC_Res_vs_NonRes = NA, P.Value = NA, adj.P.Val = NA,
      paper_direction = genes15$paper_direction[i],
      stringsAsFactors = FALSE
    )
    next
  }
  sub <- limma_res[hit, ]
  # 取 |t| 最大的探针作为该基因代表
  best <- sub[which.max(abs(sub$t)), ]
  results[[g]] <- data.frame(
    symbol = g,
    n_probe = length(hit),
    probes = paste(hit, collapse = ";"),
    best_probe = rownames(best),
    logFC_Res_vs_NonRes = round(best$logFC, 3),
    P.Value = signif(best$P.Value, 3),
    adj.P.Val = signif(best$adj.P.Val, 3),
    paper_direction = genes15$paper_direction[i]
  )
}

out <- do.call(rbind, results)

# 方向一致性判定
# GSE209746: iCR高表达 => GSE35452 Non-responder高表达 => logFC(Res vs NonRes)应为负
out$GSE35452_direction <- ifelse(out$logFC_Res_vs_NonRes > 0,
                                 "Responder(CR)高表达", "Non-responder(iCR)高表达")
out$direction_match <- with(out, {
  ours <- ifelse(logFC_Res_vs_NonRes > 0, "CR高表达", "iCR高表达")
  ifelse(is.na(logFC_Res_vs_NonRes), "无探针",
         ifelse(ours == paper_direction, "一致", "相反"))
})

cat("\n===== 15基因在 GSE35452 中的方向 =====\n")
print(out[, c("symbol", "n_probe", "best_probe", "logFC_Res_vs_NonRes",
              "P.Value", "paper_direction", "direction_match")], row.names = FALSE)

cat("\n===== 统计 =====\n")
print(table(out$direction_match))
sig_same <- sum(out$direction_match == "一致" & out$P.Value < 0.05, na.rm = TRUE)
cat("方向一致且 P<0.05 的基因数:", sig_same, "\n")

cat("\n===== IGF2 / L1CAM 全部探针 =====\n")
for (g in c("IGF2", "L1CAM")) {
  syms <- strsplit(as.character(fd[[sym_col]]), " /// ")
  hit <- fd$probe[sapply(syms, function(x) g %in% x)]
  cat("---", g, "---\n")
  print(limma_res[hit, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])
}

write.csv(out, "GSE35452_15genes_crosscheck.csv", row.names = FALSE)
cat("\n===== 完成 =====\n")
