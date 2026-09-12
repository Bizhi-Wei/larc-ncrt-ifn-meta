# ==============================
# GSE45404 队列分析 (第5队列)
# 42例 LARC nCRT, GPL570, class: Responder(19)/Non-responder(23)
# 定义: Mandard TRG (原文 PMID 26023803; PMC9121132 引述)
# 输出: features_GSE45404.csv (IFN/IGF2/L1CAM/NK/Cyto/response) + 统计
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(pROC)
})

eset <- GEOquery:::parseGSEMatrix("GSE45404-GPL570_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset)
clin <- pData(eset)
cat("样本数:", ncol(y), "\n")
cls_col <- grep("^class", colnames(clin), value = TRUE)[1]
grp <- ifelse(grepl("Non", clin[[cls_col]], ignore.case = TRUE), "poor", "good")
names(grp) <- colnames(y)
cat("分组:", paste(names(table(grp)), table(grp), collapse = " / "), "\n")

# ---- 探针→基因符号: 复用 GSE35452 (同平台 GPL570) 的注释 ----
eset_ref <- readRDS("GSE35452_eset.rds")
fd_ref <- fData(eset_ref)
ref_map <- setNames(as.character(fd_ref[["Gene symbol"]]), rownames(fd_ref))
sym_full <- ref_map[rownames(y)]
cat("探针注释覆盖率: ", round(mean(!is.na(sym_full)) * 100, 1), "%\n")
sym_first <- sapply(strsplit(sym_full, " */// *"), `[`, 1)

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）
y_sym <- symbolize(y, sym_first)

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
ifn <- colMeans(sc_ifn)

z <- function(x) as.numeric(scale(x))
get_gene <- function(g) {
  hit <- which(sapply(strsplit(sym_full, " */// *|, *"), function(x) g %in% x))
  if (length(hit) == 0) stop(g, " 无探针")
  v <- apply(y[hit, , drop = FALSE], 1, var)
  as.numeric(y[hit[which.max(v)], ])
}
feat <- data.frame(IFN = z(ifn[colnames(y)]),
                   IGF2 = z(get_gene("IGF2")),
                   L1CAM = z(get_gene("L1CAM")),
                   NK = z(sc_mcp["NK cells", colnames(y)]),
                   Cyto = z(sc_mcp["Cytotoxic lymphocytes", colnames(y)]),
                   response = grp[colnames(y)],
                   row.names = colnames(y))

# ---- 统计: Hedges g + t检验 + AUC ----
hedges <- function(x1, x2) {
  n1 <- length(x1); n2 <- length(x2)
  sp <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
  d <- (mean(x1) - mean(x2)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  se <- sqrt(J^2 * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2))))
  c(g = g, se = se)
}

cat("\n===== GSE45404 特征统计 (good=Responder) =====\n")
res <- data.frame()
for (nm in c("IFN", "IGF2", "L1CAM", "NK", "Cyto")) {
  x1 <- feat[[nm]][feat$response == "good"]
  x2 <- feat[[nm]][feat$response == "poor"]
  tt <- t.test(x1, x2)
  h <- hedges(x1, x2)
  r <- roc(feat$response, feat[[nm]], levels = c("poor", "good"),
           direction = "<", quiet = TRUE)
  cat(sprintf("%-6s g=%+.3f (se %.3f)  P=%.4f  AUC=%.3f  方向:%s\n",
              nm, h["g"], h["se"], tt$p.value, as.numeric(auc(r)),
              ifelse(mean(x1) > mean(x2), "good高", "poor高")))
  res <- rbind(res, data.frame(feature = nm, g = h["g"], se = h["se"],
                               P = tt$p.value, AUC = as.numeric(auc(r))))
}

write.csv(feat, "features_GSE45404.csv")
write.csv(res, "GSE45404_stats.csv", row.names = FALSE)

# ---- 15基因方向核对 ----
suppressMessages(library(limma))
g_f <- factor(grp, levels = c("poor", "good"))
design <- model.matrix(~ g_f)
genes15 <- data.frame(
  symbol = c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300"),
  paper_dir = c(rep("iCR_high", 8), rep("CR_high", 7)))
cat("\n===== 15基因在 GSE45404 =====\n")
n_same <- 0; n_opp <- 0; n_na <- 0
for (i in seq_len(nrow(genes15))) {
  g <- genes15$symbol[i]
  hit <- which(sapply(strsplit(sym_full, " */// *"), function(x) g %in% x))
  if (length(hit) == 0) { cat(g, ": 无探针\n"); n_na <- n_na + 1; next }
  fitp <- eBayes(lmFit(y[hit, , drop = FALSE], design))
  rp <- topTable(fitp, coef = 2, number = Inf, sort.by = "none")
  best <- rp[which.max(abs(rp$t)), ]
  d <- ifelse(best$logFC > 0, "CR_high", "iCR_high")
  same <- ifelse(d == genes15$paper_dir[i], "same", "opposite")
  if (same == "same") n_same <- n_same + 1 else n_opp <- n_opp + 1
  cat(sprintf("%-9s logFC=%+.3f P=%.3f %s(%s)\n", g, best$logFC, best$P.Value, d, same))
}
cat("一致:", n_same, " 相反:", n_opp, " 无探针:", n_na, "\n")

cat("\n===== 完成 =====\n")
