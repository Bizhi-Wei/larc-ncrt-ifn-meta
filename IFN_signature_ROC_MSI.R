# ==============================
# IFN 签名评分: 三队列 ROC/AUC
# + GSE209746 校正 MSI 后重验
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(DESeq2)
  library(pROC)
  library(limma)
})

# ---- Hallmark IFN 基因集 ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]

# 符号化表达矩阵(去重)：跨队列公共工具，见 cohort_utils.R
source("cohort_utils.R")

# ==============================
# 队列1: GSE209746
# ==============================
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym1 <- ens2sym[sub("\\..*$", "", rownames(assay(vsd)))]
y1 <- symbolize(assay(vsd), sym1)
grp1 <- ifelse(labels$response[match(colnames(y1), labels$patient)] == "CR", "good", "poor")
names(grp1) <- colnames(y1)

sc1 <- gsva(ssgseaParam(y1, ifn_sets), verbose = FALSE)
ifn1 <- colMeans(sc1)  # 两个IFN集的平均分

# ==============================
# 队列2: GSE35452
# ==============================
eset2 <- readRDS("GSE35452_eset.rds")
y2raw <- exprs(eset2)
fd2 <- fData(eset2)
clin2 <- pData(eset2)
rt <- as.character(clin2$characteristics_ch1.1)
grp2 <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
names(grp2) <- colnames(y2raw)
sym2 <- sapply(strsplit(as.character(fd2[["Gene symbol"]]), " */// *"), `[`, 1)
y2 <- symbolize(y2raw, sym2)
sc2 <- gsva(ssgseaParam(y2, ifn_sets), verbose = FALSE)
ifn2 <- colMeans(sc2)

# ==============================
# 队列3: GSE87211
# ==============================
eset3 <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
y3raw <- exprs(eset3)
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
y3t <- y3raw[, ok3]
names(grp3) <- colnames(y3t)
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym3 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y3t)]
y3 <- symbolize(y3t, sym3)
sc3 <- gsva(ssgseaParam(y3, ifn_sets), verbose = FALSE)
ifn3 <- colMeans(sc3)

# ==============================
# ROC / AUC
# ==============================
png("IFN_score_ROC_3cohorts.png", width = 1000, height = 900, res = 130)
roc_list <- list()
cohorts <- list(
  "GSE209746 (RNA-seq, n=93)"   = list(score = ifn1, mat = sc1, grp = grp1),
  "GSE35452 (Affymetrix, n=46)" = list(score = ifn2, mat = sc2, grp = grp2),
  "GSE87211 (Agilent, n=202)"   = list(score = ifn3, mat = sc3, grp = grp3)
)
cols <- c("black", "gray50", "tomato")
cat("===== IFN 签名评分 ROC/AUC =====\n")
first <- TRUE
for (i in seq_along(cohorts)) {
  nm <- names(cohorts)[i]
  sc <- cohorts[[i]]$score
  gr <- cohorts[[i]]$grp[names(sc)]
  r <- roc(gr, as.numeric(sc), levels = c("poor", "good"),
           direction = "<", quiet = TRUE)
  ci <- ci.auc(r)
  wt <- wilcox.test(sc[gr == "good"], sc[gr == "poor"])
  cat(sprintf("%-32s AUC = %.3f (95%%CI %.3f-%.3f)  Wilcoxon P = %.2e\n",
              nm, as.numeric(auc(r)), ci[1], ci[3], wt$p.value))
  plot(r, col = cols[i], lwd = 2.5, add = !first,
       main = "IFN signature score ROC", legacy.axes = TRUE)
  first <- FALSE
  roc_list[[nm]] <- r
}
abline(0, 1, lty = 3, col = "gray70")
legend("bottomright",
       legend = sapply(seq_along(cohorts), function(i)
         sprintf("%s  AUC=%.2f", names(cohorts)[i], as.numeric(auc(roc_list[[i]])))),
       col = cols, lwd = 2.5, cex = 0.85)
dev.off()

# 分开看 IFN-alpha / IFN-gamma
cat("\n===== IFN-alpha / IFN-gamma 分开的 AUC =====\n")
for (i in seq_along(cohorts)) {
  nm <- names(cohorts)[i]
  sm <- cohorts[[i]]$mat
  gr <- cohorts[[i]]$grp[colnames(sm)]
  for (s in rownames(sm)) {
    r <- roc(gr, as.numeric(sm[s, ]), levels = c("poor", "good"),
             direction = "<", quiet = TRUE)
    cat(sprintf("%-32s %-40s AUC = %.3f\n", nm, s, as.numeric(auc(r))))
  }
}

# ==============================
# GSE209746: 校正 MSI 后重验
# ==============================
cat("\n===== GSE209746 校正 MSI =====\n")
s1d <- read.csv("paper_S1D_bindea_scores.csv", stringsAsFactors = FALSE)
msi <- setNames(s1d$MSI_status, s1d$patient)[colnames(y1)]
ok_msi <- !is.na(msi)
cat("MSI 分布 (93例):", paste(names(table(msi)), table(msi), collapse = " / "),
    " 缺失:", sum(!ok_msi), "\n")
cat("MSI x response:\n")
print(table(msi[ok_msi], grp1[ok_msi]))

# 未校正 vs 校正 MSI 的 camera 对比(VST + limma, 剔除 MSI 缺失样本)
y1m <- y1[, ok_msi]
grp_f <- factor(grp1[ok_msi], levels = c("poor", "good"))
msi_f <- factor(ifelse(msi[ok_msi] == "1", "MSI", "MSS"), levels = c("MSS", "MSI"))
idx_all <- ids2indices(hallmark, rownames(y1m))

d0 <- model.matrix(~ grp_f)
d1 <- model.matrix(~ msi_f + grp_f)

cam0 <- camera(y1m, idx_all, d0, contrast = 2)
cam1 <- camera(y1m, idx_all, d1, contrast = 3)

key_sets <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_ALLOGRAFT_REJECTION",
              "HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
cat("\n通路                未校正P    校正MSI后P\n")
for (s in key_sets) {
  cat(sprintf("%-45s %.2e   %.2e\n", s, cam0[s, "PValue"], cam1[s, "PValue"]))
}

# DESeq2 校正 MSI: IGF2 / L1CAM 是否仍然显著
cat("\n===== DESeq2 ~msi+response: IGF2/L1CAM =====\n")
exp_raw <- as.matrix(read.csv(gzfile("GSE209746_geneExp_n_114_ensembl.csv.gz"),
                              row.names = 1, check.names = FALSE))
msi_lab <- setNames(s1d$MSI_status, s1d$patient)[labels$patient]
ok_lab <- !is.na(msi_lab)
lab2 <- labels[ok_lab, ]
count_mat <- exp_raw[, lab2$patient]
count_filt <- count_mat[rowSums(count_mat) >= 10, ]
coldata <- data.frame(
  response = factor(lab2$response, levels = c("CR", "iCR")),
  msi = factor(ifelse(msi_lab[ok_lab] == "1", "MSI", "MSS"), levels = c("MSS", "MSI")),
  row.names = lab2$patient
)
dds2 <- DESeqDataSetFromMatrix(count_filt, coldata, design = ~ msi + response)
dds2 <- DESeq(dds2, quiet = TRUE)
res2 <- as.data.frame(results(dds2, contrast = c("response", "iCR", "CR")))
res2$symbol <- ens2sym[sub("\\..*$", "", rownames(res2))]
chk <- res2[res2$symbol %in% c("IGF2", "L1CAM"), ]
print(chk[, c("symbol", "log2FoldChange", "pvalue", "padj")], row.names = FALSE)

cat("\n===== 完成 =====\n")
