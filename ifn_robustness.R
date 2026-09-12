# ==============================
# IFN 签名稳健性分析 (导师意见 #3)
# 问题: IFN 结果依赖某一套基因集, 还是代表稳定的生物学状态?
# 设计: primary = Hallmark IFN-α+γ ssGSEA (既有);
#       alt1 = Ayers IFN-γ 6基因 (JCI 2017, PMID 28650338);
#       alt2 = Reactome Interferon Alpha/Beta Signaling (R-HSA-909733);
#       alt3 = Reactome Interferon Gamma Signaling (R-HSA-877300);
#       key genes = CXCL9/CXCL10/IDO1/STAT1/IFNG (原始表达)
# 每队列 Hedges g (responders vs poor), REML-mKH 合并;
# 另算各替代分数与 primary IFN 分数的 Spearman 相关
# 输出: ifn_robustness_percohort.csv / ifn_robustness_pooled.csv /
#       fig_S_ifn_robustness.png/.pdf/.svg
# ==============================

suppressMessages({
  library(GEOquery); library(GSVA); library(Biobase); library(org.Hs.eg.db)
  library(illuminaHumanv3.db); library(DESeq2); library(ggplot2)
})
source("project_config.R", encoding = "UTF-8")
source("cohort_utils.R")   # symbolize(), gb2sym(), load_signature_sets()
source("meta_utils.R")     # hedges_g(), pool_random_effects()

# ---- 基因集 ----
sig <- load_signature_sets()
hallmark_ifn <- sig$ifn_sets
ayers <- list(IFNG_Ayers6 = c("IFNG", "STAT1", "IDO1", "CXCL9", "CXCL10", "HLA-DRA"))
rg <- readLines("reactome2022.gmt", warn = FALSE)
rsets <- lapply(strsplit(rg, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
rsets <- unlist(rsets, recursive = FALSE)
names(rsets) <- sapply(strsplit(rg, "\t"), `[`, 1)
react_ifn <- rsets[c("Interferon Alpha/Beta Signaling R-HSA-909733",
                     "Interferon Gamma Signaling R-HSA-877300")]
cat("Reactome IFN 集合基因数:", sapply(react_ifn, length), "\n")
alt_sets <- c(ayers, react_ifn)
key_genes <- c("CXCL9", "CXCL10", "IDO1", "STAT1", "IFNG")

# ---- 逐队列重建 symbol 矩阵 + 分组 (加载器与各原脚本逐行一致) ----
loaders <- list()

loaders$GSE209746 <- function() {
  vsd <- readRDS("GSE209746_vsd.rds")
  labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
  s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
  ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
  sym <- ens2sym[sub("\\..*$", "", rownames(assay(vsd)))]
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
  tissue_col <- grep("^tissue", colnames(clin), value = TRUE)[1]
  ypT_col <- grep("depth of invasion after", colnames(clin), value = TRUE)[1]
  ypN_col <- grep("lymph node metastasis after", colnames(clin), value = TRUE)[1]
  is_tum <- clin[[tissue_col]] == "rectal tumor"
  ypT <- suppressWarnings(as.numeric(as.character(clin[[ypT_col]])))
  ypN <- suppressWarnings(as.numeric(as.character(clin[[ypN_col]])))
  grp_all <- ifelse(ypT == 0 & ypN == 0, "good",
                    ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0), "poor", NA))
  ok <- is_tum & !is.na(grp_all)
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
  ref_map <- setNames(as.character(fData(e_ref)[["Gene symbol"]]), rownames(fData(e_ref)))
  sym <- ref_map[rownames(y)]
  list(y = symbolize(y, sym), grp = setNames(grp, colnames(y)))
}
loaders$GSE119409 <- function() {
  eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE119409_series_matrix.txt.gz",
                                    AnnotGPL = FALSE, getGPL = FALSE)$eset
  y <- exprs(eset); pd <- pData(eset)
  e_ref <- readRDS("GSE35452_eset.rds")
  ref_map <- setNames(as.character(fData(e_ref)[["Gene symbol"]]), rownames(fData(e_ref)))
  sym <- ref_map[rownames(y)]
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

ncrt9 <- LARC_CONFIG$analysis_sets$primary_ncrt9

# ---- 主循环 ----
per_rows <- list()
cor_rows <- list()
for (co in ncrt9) {
  cat("\n----", co, "----\n")
  d <- loaders[[co]]()
  y <- d$y; grp <- d$grp[colnames(y)]
  keep <- !is.na(grp)
  y <- y[, keep]; grp <- grp[keep]
  # 替代签名 ssGSEA
  sc <- gsva(ssgseaParam(y, alt_sets), verbose = FALSE)
  feats <- as.data.frame(t(sc))
  # 关键基因原始表达
  for (g in key_genes) feats[[g]] <- if (g %in% rownames(y)) as.numeric(y[g, ]) else NA_real_
  # primary IFN 分数 (既有 features 表, 按样本号匹配)
  feat0 <- read.csv(paste0("features_", co, ".csv"), row.names = 1, check.names = FALSE)
  common <- intersect(rownames(feats), rownames(feat0))
  for (ft in colnames(feats)) {
    h <- tryCatch(hedges_g(feats[[ft]], grp), error = function(e) NULL)
    if (is.null(h)) next
    per_rows[[length(per_rows) + 1]] <- data.frame(
      cohort = co, feature = ft, g = unname(h["g"]), se = unname(h["se"]))
  }
  for (ft in rownames(sc)) {
    if (length(common) >= 10)
      cor_rows[[length(cor_rows) + 1]] <- data.frame(
        cohort = co, feature = ft,
        spearman = cor(feat0[common, "IFN"], feats[common, ft], method = "spearman"))
  }
  cat(sprintf("n=%d  完成 %d 个特征\n", length(grp), ncol(feats)))
}

per <- do.call(rbind, per_rows)
cors <- do.call(rbind, cor_rows)
write.csv(per, "ifn_robustness_percohort.csv", row.names = FALSE)
write.csv(cors, "ifn_robustness_correlations.csv", row.names = FALSE)

# ---- REML-mKH 合并 ----
feats_all <- unique(per$feature)
pool_rows <- list()
for (ft in feats_all) {
  dd <- per[per$feature == ft, ]
  p <- pool_random_effects(dd$g, dd$se)
  pool_rows[[length(pool_rows) + 1]] <- data.frame(
    feature = ft, k = p$k, g = p$g, ci_lo = p$lo, ci_hi = p$hi, p = p$p,
    I2 = p$I2, tau2 = p$tau2, Q = p$Q, pi_lo = p$pi_lo, pi_hi = p$pi_hi,
    n_dir = sum(dd$g > 0))
}
pooled <- do.call(rbind, pool_rows)
# primary IFN 对照 (既有 REML-mKH 主分析结果)
primary <- read.csv("meta_pooled_results.csv", stringsAsFactors = FALSE)
p_ifn <- primary[primary$analysis == "primary_ncrt9" & primary$feature == "IFN", ]
pooled <- rbind(data.frame(feature = "IFN_Hallmark_primary", k = p_ifn$k, g = p_ifn$estimate,
                           ci_lo = p_ifn$ci_lo, ci_hi = p_ifn$ci_hi, p = p_ifn$p,
                           I2 = p_ifn$I2, tau2 = p_ifn$tau2, Q = p_ifn$Q,
                           pi_lo = p_ifn$pi_lo, pi_hi = p_ifn$pi_hi, n_dir = p_ifn$direction_consistent),
                pooled)
write.csv(pooled, "ifn_robustness_pooled.csv", row.names = FALSE)

cat("\n===== 合并结果 (9 nCRT 队列, REML-mKH) =====\n")
for (i in seq_len(nrow(pooled)))
  cat(sprintf("%-26s g=%+.2f (%.2f, %.2f) P=%.4g I2=%.0f%% 方向 %d/%d\n",
              pooled$feature[i], pooled$g[i], pooled$ci_lo[i], pooled$ci_hi[i],
              pooled$p[i], pooled$I2[i], pooled$n_dir[i], pooled$k[i]))
cat("\nSpearman 相关 (vs primary IFN):\n")
print(aggregate(spearman ~ feature, cors, function(x) round(c(min = min(x), median = median(x)), 3)))

# ---- 图: 合并效应一览 ----
lab_map <- c(IFN_Hallmark_primary = "Primary: Hallmark IFN-\u03b1/\u03b3 (ssGSEA)",
             IFNG_Ayers6 = "Alt: Ayers IFN-\u03b3 6-gene",
             `Interferon Alpha/Beta Signaling R-HSA-909733` = "Alt: Reactome IFN-\u03b1/\u03b2 signaling",
             `Interferon Gamma Signaling R-HSA-877300` = "Alt: Reactome IFN-\u03b3 signaling",
             CXCL9 = "CXCL9 (single gene)", CXCL10 = "CXCL10 (single gene)",
             IDO1 = "IDO1 (single gene)", STAT1 = "STAT1 (single gene)",
             IFNG = "IFNG (single gene)")
pooled$lab <- factor(lab_map[pooled$feature], levels = rev(unname(lab_map)))
pooled$grp_type <- ifelse(pooled$feature == "IFN_Hallmark_primary", "primary",
                          ifelse(pooled$feature %in% key_genes, "gene", "alt"))
p <- ggplot(pooled, aes(x = g, y = lab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray55") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.25) +
  geom_point(aes(color = grp_type), size = 3) +
  geom_text(aes(x = ci_hi, label = sprintf("g=%.2f (%.2f\u2013%.2f), P=%.3g",
                                           g, ci_lo, ci_hi, p)),
            hjust = -0.08, size = 3) +
  scale_color_manual(values = c(primary = "#D55E00", alt = "#0072B2", gene = "gray40"),
                     guide = "none") +
  coord_cartesian(xlim = c(min(pooled$ci_lo) - 0.1, max(pooled$ci_hi) + 1.0)) +
  labs(title = "IFN signature robustness across independent gene sets (9 nCRT cohorts)",
       x = "Pooled Hedges g (> 0: higher in responders), REML + modified Hartung-Knapp",
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        axis.text.y = element_text(color = "black"))
ggsave("fig_S_ifn_robustness.png", p, width = 10.5, height = 4.8, dpi = 300)
ggsave("fig_S_ifn_robustness.pdf", p, width = 10.5, height = 4.8,
       device = grDevices::cairo_pdf)
svglite::svglite("fig_S_ifn_robustness.svg", width = 10.5, height = 4.8)
print(p)
dev.off()
cat("\n输出: ifn_robustness_percohort.csv / ifn_robustness_pooled.csv / ",
    "ifn_robustness_correlations.csv / fig_S_ifn_robustness.png/.pdf/.svg\n",
    "===== 完成 =====\n", sep = "")
