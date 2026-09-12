# ARISTOTLE 表达预处理 + 冻结 IFN 打分 (对照臂 121)
suppressMessages({ library(oligo); library(Biobase); library(GSVA); library(clariomshumantranscriptcluster.db) })
setwd("D:/GEO2R")
source("cohort_utils.R")

cels <- list.files("scort/aristotle_cel/expression", pattern = "[.]CEL$",
                   full.names = TRUE, ignore.case = TRUE)
# 只取样本表内的对照臂 121 个
sids <- read.csv("scort/samples.csv", stringsAsFactors = FALSE)
ari <- sids$scort_id[sids$Cohort == "Aristotle"]
keep <- cels[toupper(basename(cels)) %in% paste0(ari, ".CEL")]
cat("对照臂 CEL:", length(keep), "\n")
ab <- oligo::read.celfiles(filenames = keep, pkgname = "pd.clariom.s.human", verbose = FALSE)
eset <- oligo::rma(ab, background = TRUE, normalize = TRUE)
cat("RMA:", dim(exprs(eset)), "\n")
y <- exprs(eset)
sym_vec <- unname(AnnotationDbi::mapIds(clariomshumantranscriptcluster.db,
        keys = featureNames(eset), column = "SYMBOL", keytype = "PROBEID"))
cat("symbol coverage:", round(mean(!is.na(sym_vec)) * 100, 1), "%\n")
y_sym <- symbolize(y, sym_vec)
cat("symbol matrix:", dim(y_sym), "\n")
sig <- load_signature_sets()
sc <- gsva(ssgseaParam(as.matrix(y_sym), sig$ifn_sets), verbose = FALSE)
ifn <- colMeans(sc)
out <- data.frame(sample_id = names(ifn), IFN = as.numeric(ifn))
write.csv(out, "scort/aristotle_ifn_score.csv", row.names = FALSE)
write.csv(y_sym, "scort/aristotle_expr_symbol.csv")
cat("IFN: min", round(min(out$IFN), 3), "med", round(median(out$IFN), 3),
    "max", round(max(out$IFN), 3), "\n")
