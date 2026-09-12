# ==============================
# 新队列统一处理: 表达矩阵 → symbol → IFN/NK/Cyto ssGSEA + IGF2/L1CAM
# 输出 features_GSExxxxx.csv (与既有 features 表同构)
# ==============================

suppressMessages({
  library(GEOquery); library(GSVA); library(Biobase)
})
source("cohort_utils.R")
# ---- 共享: 基因集（公共实现见 cohort_utils.R） ----
sig_sets <- load_signature_sets()
ifn_sets <- sig_sets$ifn_sets; mcp_sets <- sig_sets$mcp_sets

# symbolize() 与 build_cohort_features() 统一定义于 cohort_utils.R

# ==============================
# 1. GSE119409 (GPL570, 注释复用 GSE35452 fData)
#    response: sensitive=good / resistant=poor / unknown=drop
# ==============================
eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE119409_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset); pd <- pData(eset)
e_ref <- readRDS("GSE35452_eset.rds")
ref_map <- setNames(as.character(fData(e_ref)[["Gene symbol"]]), rownames(fData(e_ref)))
sym_full <- ref_map[rownames(y)]
sens <- sub("^[^:]*: *", "", as.character(pd[[grep("sensitivity", colnames(pd), value=TRUE)[1]]]))
grp <- ifelse(sens == "sensitive", "good", ifelse(sens == "resistant", "poor", NA))
cat("GSE119409 平台:", unique(as.character(pd$platform_id)), " 分组:", table(grp, useNA="ifany"), "\n")
y_sym <- symbolize(y, sym_full)
build_cohort_features(y_sym, setNames(grp, colnames(y)), "GSE119409", ifn_sets, mcp_sets)

# ==============================
# 2. GSE213331 (RNA-seq readcount, symbol 直接可用)
#    pCR=good / npCR=poor / Pre=drop
# ==============================
cnt <- read.delim(gzfile("new_cohorts/GSE213331_readcount_matrix.txt.gz"),
                  stringsAsFactors = FALSE, check.names = FALSE)
eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE213331_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
pd <- pData(eset)
grp_t <- setNames(sub("^group: *", "", as.character(pd[[grep("group", colnames(pd))[1]]])),
                  as.character(pd$title))
rownames(cnt) <- cnt$GeneID
cnt <- cnt[, setdiff(colnames(cnt), "GeneID")]
cnt <- cnt[, colnames(cnt) != "yangcan"]   # 无名样本剔除
cpm <- t(t(as.matrix(cnt)) / colSums(as.matrix(cnt)) * 1e6)
y <- log2(cpm + 1)
grp <- ifelse(grp_t[colnames(y)] == "pCR", "good",
              ifelse(grp_t[colnames(y)] == "npCR", "poor", NA))
names(grp) <- colnames(y)
cat("GSE213331 分组:", table(grp, useNA="ifany"), "\n")
# 去重 symbol
v <- apply(y, 1, var); o <- order(rownames(y), -v)
y <- y[o, ]; y <- y[!duplicated(rownames(y)), ]
build_cohort_features(y, grp, "GSE213331", ifn_sets, mcp_sets)

cat("===== 第一批完成 =====\n")
