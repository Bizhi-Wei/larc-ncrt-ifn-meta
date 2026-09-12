# ==============================
# ⑦ IFN 替代打分方法敏感性 (导师清单⑦)
# 正文打分 = ssGSEA(秩次)。替代 = 基因级 z-score 后对 Hallmark 两 IFN 集合
# 取均值(mean-z)。z 在队列内计算 → 仍无跨队列合并表达矩阵。
# 复用 ifn_robustness.R 已落盘的 per-cohort g? 不行——那也是 ssGSEA。
# 本脚本直接重建 9 队列 symbol 矩阵, mean-z 打分, hedges_g, REML-mKH。
# 输出: ifn_scoring_sensitivity.csv / ifn_scoring_pooled.csv
# ==============================

suppressMessages({
  library(GEOquery); library(Biobase); library(org.Hs.eg.db)
  library(illuminaHumanv3.db); library(DESeq2)
})
setwd("D:/GEO2R")
source("cohort_utils.R")
source("meta_utils.R")

sig <- load_signature_sets()
ifn_genes <- unique(unlist(sig$ifn_sets))
cat("Hallmark IFN 两集合合并基因数:", length(ifn_genes), "\n")

# ---- 与 mcp_package_sensitivity.R 相同的 9 队列加载器 ----
loaders <- list()
loaders$GSE209746 <- function() {
  vsd <- readRDS("GSE209746_vsd.rds")
  labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
  s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
  sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)[sub("\\..*$", "", rownames(assay(vsd)))]
  y <- symbolize(assay(vsd), sym)
  grp <- ifelse(labels$response[match(colnames(y), labels$patient)] == "CR", "good", "poor")
  list(y = y, grp = setNames(grp, colnames(y)))
}
loaders$GSE35452 <- function() {
  eset <- readRDS("GSE35452_eset.rds")
  y <- exprs(eset); clin <- pData(eset); fd <- fData(eset)
  rt <- as.character(clin$characteristics_ch1.1)
  grp <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
  sym <- sapply(strsplit(as.character(fd[["Gene symbol"]]), " */// *"), `[`, 1)
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE87211 <- function() {
  eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); clin <- pData(eset)
  tis <- clin[[grep("^tissue", colnames(clin), value = TRUE)[1]]]
  ypT <- suppressWarnings(as.numeric(as.character(clin[[grep("depth of invasion after", colnames(clin), value = TRUE)[1]]])))
  ypN <- suppressWarnings(as.numeric(as.character(clin[[grep("lymph node metastasis after", colnames(clin), value = TRUE)[1]]])))
  grp_all <- ifelse(ypT == 0 & ypN == 0, "good",
                    ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0), "poor", NA))
  ok <- tis == "rectal tumor" & !is.na(grp_all)
  grp <- grp_all[ok]
  stopifnot(sum(grp == "good") == 35L, sum(grp == "poor") == 167L)
  yt <- y[, ok]
  pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
  sym <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(yt)]
  list(y = symbolize(yt, sym), grp = setNames(grp, colnames(yt)))
}
loaders$GSE150082 <- function() {
  eset <- GEOquery:::parseGSEMatrix("GSE150082_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); clin <- pData(eset)
  grp <- ifelse(clin$`response:ch1` == "Good", "good", "poor")
  pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
  sym <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y)]
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE45404 <- function() {
  eset <- GEOquery:::parseGSEMatrix("GSE45404-GPL570_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); clin <- pData(eset)
  cls_col <- grep("^class", colnames(clin), value = TRUE)[1]
  grp <- ifelse(grepl("Non", clin[[cls_col]], ignore.case = TRUE), "poor", "good")
  e_ref <- readRDS("GSE35452_eset.rds")
  sym <- setNames(as.character(fData(e_ref)[["Gene symbol"]]), rownames(fData(e_ref)))[rownames(y)]
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE119409 <- function() {
  eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE119409_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); pd <- pData(eset)
  e_ref <- readRDS("GSE35452_eset.rds")
  sym <- setNames(as.character(fData(e_ref)[["Gene symbol"]]), rownames(fData(e_ref)))[rownames(y)]
  sens <- sub("^[^:]*: *", "", as.character(pd[[grep("sensitivity", colnames(pd), value = TRUE)[1]]]))
  grp <- ifelse(sens == "sensitive", "good", ifelse(sens == "resistant", "poor", NA))
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE133057 <- function() {
  eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE133057_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); pd <- pData(eset)
  ann <- read.delim(gzfile("new_cohorts/GPL6102.annot.gz"), skip = 28,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
  sym <- setNames(ann[["Gene.symbol"]], ann$ID)[rownames(y)]
  gbacc <- setNames(ann[["GenBank.Accession"]], ann$ID)[rownames(y)]
  need <- is.na(sym) | sym == ""
  sym[need] <- gb2sym(gbacc[need])
  ajcc <- as.numeric(sub("^[^:]*: *", "", as.character(pd[[grep("ajcc", colnames(pd), value = TRUE)[1]]])))
  grp <- ifelse(ajcc %in% c(0, 1), "good", ifelse(ajcc %in% c(2, 3), "poor", NA))
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE53781 <- function() {
  eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE53781_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); pd <- pData(eset)
  lines <- readLines(gzfile("new_cohorts/GPL18134.soft.gz"))
  b <- grep("^!platform_table_begin", lines); e <- grep("^!platform_table_end", lines)
  plat <- read.delim(textConnection(lines[(b + 1L):(e - 1L)]),
                     stringsAsFactors = FALSE, quote = "")
  sym <- unname(gb2sym(setNames(plat$GB_ACC, plat$ID)[rownames(y)]))
  resp_col <- grep("^resp", colnames(pd), value = TRUE)[1]
  resp <- sub("^[^:]*: *", "", as.character(pd[[resp_col]]))
  grp <- ifelse(resp == "YES", "good", ifelse(resp == "NO", "poor", NA))
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE94104 <- function() {
  eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE94104_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); pd <- pData(eset)
  m <- suppressMessages(AnnotationDbi::select(illuminaHumanv3.db,
          keys = rownames(y), columns = "SYMBOL", keytype = "PROBEID"))
  m <- m[!duplicated(m$PROBEID), ]
  sym <- unname(setNames(m$SYMBOL, m$PROBEID)[rownames(y)])
  ttype <- sub("^[^:]*: *", "", as.character(pd[[grep("tissue type", colnames(pd), value = TRUE)[1]]]))
  trg <- as.numeric(sub("^[^:]*: *", "", as.character(pd[[grep("tumour regression", colnames(pd), value = TRUE)[1]]])))
  is_pre <- ttype == "Biopsy"
  grp <- ifelse(trg == 3, "good", ifelse(trg %in% c(1, 2), "poor", NA))
  stopifnot(sum(grp[is_pre] == "good", na.rm = TRUE) == 11L,
            sum(grp[is_pre] == "poor", na.rm = TRUE) == 29L)
  yt <- y[, is_pre]
  list(y = symbolize(yt, sym), grp = setNames(grp[is_pre], colnames(yt)))
}

ncrt9 <- c("GSE209746", "GSE35452", "GSE87211", "GSE150082", "GSE45404",
           "GSE119409", "GSE133057", "GSE53781", "GSE94104")

rows <- list()
for (co in ncrt9) {
  d <- loaders[[co]]()
  y <- d$y; grp <- d$grp[colnames(y)]
  keep <- !is.na(grp); y <- y[, keep]; grp <- grp[keep]
  hit <- intersect(ifn_genes, rownames(y))
  # 队列内 z-score (基因级), 对缺失基因静默跳过
  yz <- y[hit, , drop = FALSE]
  yz <- t(scale(t(yz)))
  score <- colMeans(yz, na.rm = TRUE)
  h <- tryCatch(hedges_g(as.numeric(score), grp), error = function(e) NULL)
  if (is.null(h)) { cat(co, "FAIL\n"); next }
  rows[[length(rows) + 1]] <- data.frame(cohort = co, method = "mean_z",
                                         n_genes = length(hit),
                                         g = unname(h["g"]), se = unname(h["se"]))
  cat(sprintf("%-10s genes=%3d/%d  g=%+.2f\n", co, length(hit), length(ifn_genes), h["g"]))
  # 与 ssGSEA 版相关
  f0 <- read.csv(paste0("features_", co, ".csv"), row.names = 1)
  common <- intersect(names(score), rownames(f0))
  if (length(common) >= 10) {
    r <- cor(score[common], f0[common, "IFN"], method = "spearman")
    rows[[length(rows) + 1]] <- data.frame(cohort = co, method = "cor_vs_ssGSEA",
                                           n_genes = length(hit), g = r, se = NA)
  }
}
per <- do.call(rbind, rows)
write.csv(per, "ifn_scoring_sensitivity.csv", row.names = FALSE)

dd <- per[per$method == "mean_z", ]
p <- pool_random_effects(dd$g, dd$se)
pooled <- data.frame(method = "mean_z", k = p$k, g = p$g, ci_lo = p$lo,
                     ci_hi = p$hi, p = p$p, I2 = p$I2, tau2 = p$tau2,
                     pi_lo = p$pi_lo, pi_hi = p$pi_hi, n_dir = sum(dd$g > 0))
primary <- read.csv("meta_pooled_results.csv", stringsAsFactors = FALSE)
r0 <- primary[primary$analysis == "primary_ncrt9" & primary$feature == "IFN", ]
pooled <- rbind(pooled, data.frame(method = "ssGSEA [primary ref]", k = r0$k,
  g = r0$estimate, ci_lo = r0$ci_lo, ci_hi = r0$ci_hi, p = r0$p, I2 = r0$I2,
  tau2 = r0$tau2, pi_lo = r0$pi_lo, pi_hi = r0$pi_hi,
  n_dir = r0$direction_consistent))
write.csv(pooled, "ifn_scoring_pooled.csv", row.names = FALSE)

cat("\n===== 打分方法敏感性 (9 nCRT 队列) =====\n")
for (i in seq_len(nrow(pooled)))
  cat(sprintf("%-22s g=%+.2f (%.2f, %.2f) P=%.4g I2=%.0f%% 方向 %d/%d\n",
              pooled$method[i], pooled$g[i], pooled$ci_lo[i], pooled$ci_hi[i],
              pooled$p[i], pooled$I2[i], pooled$n_dir[i], pooled$k[i]))
cat("median spearman(mean-z, ssGSEA) =",
    median(per$g[per$method == "cor_vs_ssGSEA"]), "\n")
cat("输出: ifn_scoring_sensitivity.csv / ifn_scoring_pooled.csv\n===== 完成 =====\n")
