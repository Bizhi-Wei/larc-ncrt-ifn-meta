# ==============================
# 第三批: GSE53781 + GSE94104 (illuminaHumanv3.db 注释) +
#          GSE133057 注释补强重跑
# ==============================

suppressMessages({
  library(GEOquery); library(GSVA); library(Biobase); library(org.Hs.eg.db)
  library(illuminaHumanv3.db)
})
source("cohort_utils.R")
sig_sets <- load_signature_sets()
ifn_sets <- sig_sets$ifn_sets; mcp_sets <- sig_sets$mcp_sets

# symbolize() / gb2sym() / build_cohort_features() 统一定义于 cohort_utils.R

# ==============================
# 1. GSE53781 (CodeLink; GB_ACC -> symbol)
# ==============================
eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE53781_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset); pd <- pData(eset)
lines <- readLines(gzfile("new_cohorts/GPL18134.soft.gz"))
b <- grep("^!platform_table_begin", lines); e <- grep("^!platform_table_end", lines)
stopifnot(length(b) == 1L, length(e) == 1L, e > b)
plat <- read.delim(textConnection(lines[(b + 1L):(e - 1L)]),
                   stringsAsFactors = FALSE, quote = "")
gb <- setNames(plat$GB_ACC, plat$ID)[rownames(y)]
sym_full <- unname(gb2sym(gb))
cat("GSE53781 注释覆盖:", mean(!is.na(sym_full) & sym_full != ""), "\n")
resp_col <- grep("^resp", colnames(pd), value = TRUE)[1]
stopifnot(length(resp_col) == 1L)
resp <- sub("^[^:]*: *", "", as.character(pd[[resp_col]]))
grp <- ifelse(resp == "YES", "good", ifelse(resp == "NO", "poor", NA))
cat("GSE53781 分组:", paste(table(grp, useNA = "ifany"), collapse = "/"), "\n")
y_sym <- symbolize(y, sym_full)
build_cohort_features(y_sym, setNames(grp, colnames(y)), "GSE53781", ifn_sets, mcp_sets)

# ==============================
# 2. GSE94104 (Illumina HT-12 v3; 仅 pre-therapeutic biopsy)
# 原文使用升序退缩分级；GEO 活检中 TRG3 为最佳可用反应组，TRG1/2 为较差反应组。
# ==============================
eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE94104_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset); pd <- pData(eset)
m <- suppressMessages(AnnotationDbi::select(illuminaHumanv3.db,
        keys = rownames(y), columns = "SYMBOL", keytype = "PROBEID"))
m <- m[!duplicated(m$PROBEID), ]
sym_full <- unname(setNames(m$SYMBOL, m$PROBEID)[rownames(y)])
cat("GSE94104 注释覆盖:", mean(!is.na(sym_full)), "\n")
ttype <- sub("^[^:]*: *", "", as.character(pd[[grep("tissue type", colnames(pd), value=TRUE)[1]]]))
trg   <- as.numeric(sub("^[^:]*: *", "", as.character(pd[[grep("tumour regression", colnames(pd), value=TRUE)[1]]])))
cat("GSE94104 样本类型:", paste(table(ttype), collapse="/"), " TRG分布:", paste(table(trg, useNA="ifany"), collapse="/"), "\n")
is_pre <- ttype == "Biopsy"
grp <- ifelse(trg == 3, "good", ifelse(trg %in% c(1, 2), "poor", NA))
cat("活检中 TRG3:", sum(grp[is_pre] == "good", na.rm=TRUE),
    " TRG1/2:", sum(grp[is_pre] == "poor", na.rm=TRUE), "\n")
stopifnot(sum(grp[is_pre] == "good", na.rm = TRUE) == 11L,
          sum(grp[is_pre] == "poor", na.rm = TRUE) == 29L)
y_sym <- symbolize(y[, is_pre], sym_full)
build_cohort_features(y_sym, setNames(grp[is_pre], colnames(y[, is_pre])), "GSE94104", ifn_sets, mcp_sets)

# ==============================
# 3. GSE133057 注释补强 (GPL6102 annot 的 Gene symbol + GenBank 回补)
# ==============================
eset <- GEOquery:::parseGSEMatrix("new_cohorts/GSE133057_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset); pd <- pData(eset)
ann <- read.delim(gzfile("new_cohorts/GPL6102.annot.gz"), skip = 28,
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
sym0 <- setNames(ann[["Gene.symbol"]], ann$ID)[rownames(y)]
gbacc <- setNames(ann[["GenBank.Accession"]], ann$ID)[rownames(y)]
need <- is.na(sym0) | sym0 == ""
sym0[need] <- gb2sym(gbacc[need])
cat("GSE133057 回补后注释覆盖:", mean(!is.na(sym0) & sym0 != ""), "\n")
ajcc <- as.numeric(sub("^[^:]*: *", "", as.character(pd[[grep("ajcc", colnames(pd), value=TRUE)[1]]])))
grp <- ifelse(ajcc %in% c(0, 1), "good", ifelse(ajcc %in% c(2, 3), "poor", NA))
y_sym <- symbolize(y, sym0)
build_cohort_features(y_sym, setNames(grp, colnames(y)), "GSE133057", ifn_sets, mcp_sets)

cat("===== 第三批完成 =====\n")
