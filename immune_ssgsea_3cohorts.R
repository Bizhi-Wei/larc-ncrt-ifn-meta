# ==============================
# 三队列免疫细胞浸润 ssGSEA 分析
# MCP-counter 10类细胞签名, 方向统一: good = CR/pCR/Responder
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
})

# ---- MCP-counter 签名 ----
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
cat("MCP-counter 细胞类型数:", length(mcp_sets), "\n")
print(sapply(mcp_sets, length))

# ---- 通用: ssGSEA + 组间 t 检验 ----
run_ssgsea <- function(y, symbols, group) {
  keep <- !is.na(symbols) & symbols != ""
  y <- y[keep, ]; symbols <- symbols[keep]
  v <- apply(y, 1, var)
  ord <- order(symbols, -v)
  y <- y[ord, ]; symbols <- symbols[ord]
  dup <- duplicated(symbols)
  y <- y[!dup, ]; rownames(y) <- symbols[!dup]

  par <- ssgseaParam(as.matrix(y), mcp_sets)
  sc <- gsva(par, verbose = FALSE)

  res <- data.frame()
  for (ct in rownames(sc)) {
    g1 <- sc[ct, group == "good"]
    g2 <- sc[ct, group == "poor"]
    tt <- t.test(g1, g2)
    res <- rbind(res, data.frame(
      cell = ct,
      good_mean = round(mean(g1), 4),
      poor_mean = round(mean(g2), 4),
      direction = ifelse(mean(g1) > mean(g2), "Up_in_good", "Up_in_poor"),
      P = tt$p.value
    ))
  }
  res[order(res$P), ]
}

# ---- 队列1: GSE209746 (VST) ----
suppressMessages(library(DESeq2))
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
y1 <- assay(vsd)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym1 <- ens2sym[sub("\\..*$", "", rownames(y1))]
grp1 <- ifelse(labels$response[match(colnames(y1), labels$patient)] == "CR", "good", "poor")
r1 <- run_ssgsea(y1, sym1, grp1)
colnames(r1)[2:5] <- paste0(colnames(r1)[2:5], "_GSE209746")

# ---- 队列2: GSE35452 ----
eset2 <- readRDS("GSE35452_eset.rds")
y2 <- exprs(eset2)
fd2 <- fData(eset2)
clin2 <- pData(eset2)
rt <- as.character(clin2$characteristics_ch1.1)
grp2 <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
sym2 <- sapply(strsplit(as.character(fd2[["Gene symbol"]]), " */// *"), `[`, 1)
r2 <- run_ssgsea(y2, sym2, grp2)
colnames(r2)[2:5] <- paste0(colnames(r2)[2:5], "_GSE35452")

# ---- 队列3: GSE87211 (仅肿瘤) ----
eset3 <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
y3 <- exprs(eset3)
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
y3 <- y3[, ok3]
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym3 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y3)]
r3 <- run_ssgsea(y3, sym3, grp3)
colnames(r3)[2:5] <- paste0(colnames(r3)[2:5], "_GSE87211")

# ---- 合并 ----
m <- merge(merge(r1, r2, by = "cell"), r3, by = "cell")
m$consistent <- (m$direction_GSE209746 == m$direction_GSE35452) &
                (m$direction_GSE35452 == m$direction_GSE87211)
m$min_P <- pmin(m$P_GSE209746, m$P_GSE35452, m$P_GSE87211)
m <- m[order(m$P_GSE209746 * m$P_GSE35452 * m$P_GSE87211), ]

cat("\n===== 三队列免疫细胞 ssGSEA 对照 =====\n")
out <- m[, c("cell",
             "direction_GSE209746", "P_GSE209746",
             "direction_GSE35452", "P_GSE35452",
             "direction_GSE87211", "P_GSE87211", "consistent")]
for (cc in grep("^P_", colnames(out), value = TRUE)) out[[cc]] <- signif(out[[cc]], 3)
print(out, row.names = FALSE)
write.csv(m, "immune_ssgsea_3cohorts.csv", row.names = FALSE)

# ---- GSE209746: 我们的 MCP 分数 vs 论文 BINDEA 分数(方法学验证) ----
vsd2 <- readRDS("GSE209746_vsd.rds")
par1 <- ssgseaParam(as.matrix(y1_dup <- {  # 复用 run_ssgsea 内部处理的矩阵
  keep <- !is.na(sym1) & sym1 != ""
  yy <- y1[keep, ]; ss <- sym1[keep]
  v <- apply(yy, 1, var); o <- order(ss, -v)
  yy <- yy[o, ]; ss <- ss[o]; dd <- duplicated(ss)
  rownames(yy) <- ss[!dd]; yy[!dd, ]
}), mcp_sets)
sc1 <- gsva(par1, verbose = FALSE)

bindea <- read.csv("paper_S1D_bindea_scores.csv", stringsAsFactors = FALSE)
bindea <- bindea[match(colnames(sc1), bindea$patient), ]
pairs <- list(
  c("T cells", "T_CELLS"), c("CD8 T cells", "CD8_T_CELLS"),
  c("Cytotoxic lymphocytes", "CYTOTOXIC_CELLS"), c("NK cells", "NK_CELLS"),
  c("B lineage", "B_CELLS"), c("Myeloid dendritic cells", "DC"),
  c("Neutrophils", "NEUTROPHILS")
)
cat("\n===== 方法学验证: 我们的MCP ssGSEA vs 论文BINDEA (Pearson相关) =====\n")
for (p in pairs) {
  ours <- sc1[p[1], ]
  theirs <- as.numeric(bindea[[p[2]]])
  cc <- cor(ours, theirs, use = "complete.obs")
  cat(sprintf("%-24s vs %-16s r = %.3f\n", p[1], p[2], cc))
}

cat("\n===== 完成 =====\n")
