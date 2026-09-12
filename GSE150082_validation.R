# ==============================
# GSE150082 验证: 39例 LARC 治疗前活检
# response: Good(16) / Poor(23), 平台 GPL13497(同GSE87211)
# 内容: Hallmark camera / 15基因方向 / MCP-counter / IFN AUC / 联合模型迁移
# ==============================

suppressMessages({
  library(GEOquery)
  library(limma)
  library(GSVA)
  library(pROC)
})

eset <- GEOquery:::parseGSEMatrix("GSE150082_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset)
clin <- pData(eset)
grp <- ifelse(clin$`response:ch1` == "Good", "good", "poor")
names(grp) <- colnames(y)
cat("分组:", paste(names(table(grp)), table(grp), collapse = " / "), "\n")

pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym_full <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y)]
sym_first <- sapply(strsplit(sym_full, " */// *|, *"), `[`, 1)

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）
y_sym <- symbolize(y, sym_first)

# ---- 1. Hallmark camera ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])

g_f <- factor(grp, levels = c("poor", "good"))
design <- model.matrix(~ g_f)
idx <- ids2indices(hallmark, rownames(y_sym))
cam <- camera(y_sym, idx, design, contrast = 2)

key_sets <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_ALLOGRAFT_REJECTION",
              "HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
              "HALLMARK_GLYCOLYSIS",
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_HYPOXIA")
cat("\n===== 关注通路在 GSE150082 =====\n")
for (s in key_sets) {
  cat(sprintf("%-45s %-4s P=%.2e FDR=%.2e\n",
              s, cam[s, "Direction"], cam[s, "PValue"], cam[s, "FDR"]))
}
cam$set <- rownames(cam)
write.csv(cam, "GSE150082_hallmark_camera.csv", row.names = FALSE)

# ---- 2. 15基因方向 ----
fit <- eBayes(lmFit(y_sym, design))
res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

genes15 <- data.frame(
  symbol = c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300"),
  paper_dir = c(rep("iCR_high", 8), rep("CR_high", 7)))
cat("\n===== 15基因在 GSE150082 =====\n")
n_same <- 0; n_opp <- 0; n_na <- 0
for (i in seq_len(nrow(genes15))) {
  g <- genes15$symbol[i]
  hit <- which(sapply(strsplit(sym_full, " */// *|, *"), function(x) g %in% x))
  if (length(hit) == 0) { cat(g, ": 无探针\n"); n_na <- n_na + 1; next }
  # 用原始探针级 limma
  fitp <- eBayes(lmFit(y[hit, , drop = FALSE], design))
  rp <- topTable(fitp, coef = 2, number = Inf, sort.by = "none")
  best <- rp[which.max(abs(rp$t)), ]
  d <- ifelse(best$logFC > 0, "CR_high", "iCR_high")
  same <- ifelse(d == genes15$paper_dir[i], "same", "opposite")
  if (same == "same") n_same <- n_same + 1 else n_opp <- n_opp + 1
  cat(sprintf("%-9s logFC=%+.3f P=%.3f %s(%s)\n",
              g, best$logFC, best$P.Value, d, same))
}
cat("一致:", n_same, " 相反:", n_opp, " 无探针:", n_na, "\n")

# ---- 3. MCP-counter 免疫细胞 ----
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
sc_mcp <- gsva(ssgseaParam(y_sym, mcp_sets), verbose = FALSE)
cat("\n===== MCP-counter 免疫细胞 =====\n")
for (ct in rownames(sc_mcp)) {
  g1 <- sc_mcp[ct, grp == "good"]; g2 <- sc_mcp[ct, grp == "poor"]
  tt <- t.test(g1, g2)
  cat(sprintf("%-24s %s P=%.3f\n", ct,
              ifelse(mean(g1) > mean(g2), "Up_in_good ", "Up_in_poor "), tt$p.value))
}

# ---- 4. IFN 评分 AUC + 联合模型迁移 ----
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
ifn <- colMeans(sc_ifn)

z <- function(x) as.numeric(scale(x))
get_gene <- function(g) {
  hit <- which(sapply(strsplit(sym_full, " */// *|, *"), function(x) g %in% x))
  v <- apply(y[hit, , drop = FALSE], 1, var)
  as.numeric(y[hit[which.max(v)], ])
}
feat <- data.frame(IFN = z(ifn[colnames(y)]),
                   IGF2 = z(get_gene("IGF2")),
                   L1CAM = z(get_gene("L1CAM")),
                   NK = z(sc_mcp["NK cells", ]),
                   Cyto = z(sc_mcp["Cytotoxic lymphocytes", ]),
                   response = grp[colnames(y)],
                   row.names = colnames(y))

r <- roc(feat$response, feat$IFN, levels = c("poor", "good"),
         direction = "<", quiet = TRUE)
ci <- ci.auc(r)
cat(sprintf("\nIFN 评分 AUC(GSE150082): %.3f (%.3f-%.3f)\n",
            as.numeric(auc(r)), ci[1], ci[3]))

# 联合模型系数(GSE209746 训练的固定系数, 截距无关AUC)
lp <- -0.2988 * feat$IFN + 0.0402 * feat$IGF2 + 0.5866 * feat$L1CAM
# 注意: 模型预测 poor=1, 方向取反评估
r2 <- roc(feat$response, -lp, levels = c("poor", "good"),
          direction = "<", quiet = TRUE)
cat(sprintf("联合模型(IFN+IGF2+L1CAM)迁移 AUC: %.3f\n", as.numeric(auc(r2))))

write.csv(feat, "features_GSE150082.csv")
cat("\n===== 完成 =====\n")
