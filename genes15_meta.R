# ==============================
# 15 基因逐基因跨队列 meta 分析 (Hedges g, REML-mKH 随机效应)
# g 编码: 反应组(good) vs 差反应组(poor), 正值=反应组高表达
# ==============================

suppressMessages({ library(Biobase); library(GEOquery); library(DESeq2)
                   library(ggplot2) })
source("meta_utils.R")

genes15 <- c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300")
# 论文方向(iCR_high=耐药高表达, 期望 g<0; CR_high 期望 g>0)
paper_dir <- c(VAX1="iCR_high", REG3A="iCR_high", IGF2="iCR_high",
               NHLH2="iCR_high", SOX2="iCR_high", KLK6="iCR_high",
               L1CAM="iCR_high", INHBB="iCR_high", CPS1="CR_high",
               ALOX12B="iCR_high", TDRD1="CR_high", UNC13A="iCR_high",
               TENM1="CR_high", ADCY2="iCR_high", ZNF300="CR_high")

# 取基因: 多探针取方差最大; symbol 列支持 "///" 分隔的别名
get_gene_expr <- function(expr, sym_full, gene) {
  hit <- which(sapply(strsplit(sym_full, " */// *|, *"), function(x) gene %in% x))
  if (length(hit) == 0) return(NULL)
  v <- apply(expr[hit, , drop = FALSE], 1, var)
  as.numeric(expr[hit[which.max(v)], ])
}

per_cohort <- list()

# ---- GSE209746 (RNA-seq vsd, 行名=Ensembl, 列名=JM 患者号) ----
vsd <- readRDS("GSE209746_vsd.rds")
lab <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8 <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8$hgnc_symbol, sub("\\..*", "", s8$ensembl_gene_id))
ex <- assay(vsd)
sym209 <- ens2sym[sub("\\..*", "", rownames(ex))]
keep <- !is.na(sym209)
ex <- ex[keep, ]; sym209 <- sym209[keep]
v <- apply(ex, 1, var); o <- order(sym209, -v)
ex <- ex[o, ]; sym209 <- sym209[o]
ex <- ex[!duplicated(sym209), ]
sym209 <- sym209[!duplicated(sym209)]
grp <- ifelse(lab$response[match(colnames(ex), lab$patient)] == "CR", "good", "poor")
per_cohort[["GSE209746"]] <- list(expr = ex, symbols = sym209, group = grp)

# ---- GSE35452 (GPL570 eset) ----
e <- readRDS("GSE35452_eset.rds")
ex <- exprs(e)
sym_full <- as.character(fData(e)[["Gene symbol"]])[match(rownames(ex), rownames(fData(e)))]
pd <- pData(e)
rsp_col <- grep("response|responder", colnames(pd), value = TRUE, ignore.case = TRUE)[1]
grp_raw <- as.character(pd[[rsp_col]])
grp <- ifelse(grepl("^R|responder", grp_raw, ignore.case = TRUE) &
              !grepl("non", grp_raw, ignore.case = TRUE), "good", "poor")
cat("GSE35452 response 字段:", rsp_col, " 取值:", paste(unique(grp_raw), collapse = "/"),
    " -> good:", sum(grp == "good"), "\n")
per_cohort[["GSE35452"]] <- list(expr = ex, symbols = sym_full, group = grp)

# ---- GSE87211 (GPL13497 单色) ----
eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
ex <- exprs(eset); clin <- pData(eset)
is_tum <- clin[["tissue:ch1"]] == "rectal tumor"
ypT <- suppressWarnings(as.numeric(as.character(clin[["depth of invasion after rct:ch1"]])))
ypN <- suppressWarnings(as.numeric(as.character(clin[["lymph node metastasis after rct:ch1"]])))
grp <- ifelse(ypT == 0 & ypN == 0, "good",
       ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0), "poor", NA))
ok <- is_tum & !is.na(grp)
stopifnot(sum(grp[ok] == "good") == 35L, sum(grp[ok] == "poor") == 167L)
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym_full <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(ex)]
per_cohort[["GSE87211"]] <- list(expr = ex[, ok], symbols = sym_full, group = grp[ok])

# ---- GSE150082 (GPL13497 双色) ----
eset <- GEOquery:::parseGSEMatrix("GSE150082_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
ex <- exprs(eset); clin <- pData(eset)
rsp_col <- grep("response|trg|regression", colnames(clin), value = TRUE, ignore.case = TRUE)
grp_raw <- apply(clin[, rsp_col, drop = FALSE], 1, paste, collapse = " ")
grp <- ifelse(grepl("good", grp_raw, ignore.case = TRUE), "good",
       ifelse(grepl("poor", grp_raw, ignore.case = TRUE), "poor", NA))
cat("GSE150082 response:", paste(table(grp, useNA = "ifany"), collapse = "/"), "\n")
ok <- !is.na(grp)
sym_full <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(ex)]
per_cohort[["GSE150082"]] <- list(expr = ex[, ok], symbols = sym_full, group = grp[ok])

# ---- GSE45404 (GPL570, 复用 GSE35452 注释) ----
eset <- GEOquery:::parseGSEMatrix("GSE45404-GPL570_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
ex <- exprs(eset); clin <- pData(eset)
e_ref <- readRDS("GSE35452_eset.rds")
fd_ref <- fData(e_ref)
ref_map <- setNames(as.character(fd_ref[["Gene symbol"]]), rownames(fd_ref))
sym_full <- ref_map[rownames(ex)]
rsp_col <- grep("response|trg|mandard", colnames(clin), value = TRUE, ignore.case = TRUE)
grp_raw <- apply(clin[, rsp_col, drop = FALSE], 1, paste, collapse = " ")
grp <- ifelse(grepl("responder", grp_raw, ignore.case = TRUE) &
              !grepl("non", grp_raw, ignore.case = TRUE), "good",
       ifelse(grepl("non", grp_raw, ignore.case = TRUE), "poor", NA))
if (all(is.na(grp))) {  # 退回已知分组: 前19 R / 后23 NR 不可靠, 直接读既有特征表对齐
  f <- read.csv("features_GSE45404.csv", stringsAsFactors = FALSE)
  grp <- setNames(f$response, f[[1]])[colnames(ex)]
}
cat("GSE45404 response:", paste(table(grp, useNA = "ifany"), collapse = "/"), "\n")
ok <- !is.na(grp)
sym_full[is.na(sym_full)] <- ""
per_cohort[["GSE45404"]] <- list(expr = ex[, ok], symbols = sym_full, group = grp[ok])

# ---- 逐基因逐队列 g ----
long <- list()
for (co in names(per_cohort)) {
  pc <- per_cohort[[co]]
  for (g in genes15) {
    x <- get_gene_expr(pc$expr, pc$symbols, g)
    if (is.null(x)) {
      long[[length(long) + 1]] <- data.frame(cohort = co, gene = g, g = NA, se = NA)
      next
    }
    h <- hedges_g(x, pc$group)
    long[[length(long) + 1]] <- data.frame(cohort = co, gene = g, g = h["g"], se = h["se"])
  }
}
long <- do.call(rbind, long)
write.csv(long, "genes15_meta_percohort.csv", row.names = FALSE)

# ---- REML + modified Hartung-Knapp 合并 ----
summ <- list()
for (g in genes15) {
  dd <- long[long$gene == g & !is.na(long$g), ]
  if (nrow(dd) < 3) next
  p <- pool_random_effects(dd$g, dd$se)
  dir_exp <- ifelse(paper_dir[g] == "iCR_high", -1, 1)
  summ[[g]] <- data.frame(gene = g, n_cohorts = nrow(dd), paper_dir = paper_dir[g],
                          pooled_g = p$g, lo = p$lo, hi = p$hi, p = p$p, I2 = p$I2,
                          tau2 = p$tau2, pi_lo = p$pi_lo, pi_hi = p$pi_hi,
                          method = p$method,
                          direction_match = sign(p$g) == dir_exp)
}
summ <- do.call(rbind, summ)
summ$fdr <- p.adjust(summ$p, method = "BH")
summ <- summ[order(summ$pooled_g), ]
write.csv(summ, "genes15_meta_pooled.csv", row.names = FALSE)

cat("\n===== 15基因 meta 合并 =====\n")
print(summ[, c("gene", "n_cohorts", "paper_dir", "pooled_g", "p", "fdr", "I2", "direction_match")],
      row.names = FALSE, digits = 2)
cat("方向一致:", sum(summ$direction_match), "/", nrow(summ),
    " FDR<0.05:", sum(summ$fdr < 0.05), "\n")

# ---- 森林图(合并 g, 按效应排序) ----
summ$gene <- factor(summ$gene, levels = summ$gene)
summ$match_col <- ifelse(summ$direction_match, "direction consistent", "opposite")
p_gene <- ggplot(summ, aes(x = pooled_g, y = gene)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, color = "grey60") +
  geom_point(aes(color = match_col), size = 2.6) +
  scale_color_manual(values = c("direction consistent" = "#0072B2", "opposite" = "#D55E00"),
                     name = NULL) +
  labs(title = "15-gene signature: cross-cohort random-effects meta-analysis",
       x = "Pooled Hedges g (responders vs poor responders)", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom",
                                   plot.title = element_text(size = 11, face = "bold"))
ggsave("fig_S_genes15_meta.png", p_gene, width = 6.5, height = 5.5, dpi = 300)
ggsave("fig_S_genes15_meta.pdf", p_gene, width = 6.5, height = 5.5,
       device = grDevices::cairo_pdf)

cat("\n输出: genes15_meta_percohort.csv / genes15_meta_pooled.csv / fig_S_genes15_meta.png/.pdf\n===== 完成 =====\n")
