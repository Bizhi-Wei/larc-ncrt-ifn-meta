# ==============================
# ⑤-a Grampian 表达预处理 + 冻结 IFN 打分
# 输入: scort/grampian.zip (223 Affymetrix Xcel/CEL, Affy 官方样本表全匹配)
# 输出: scort/grampian_expr_symbol.csv (log2 RMA, symbol 行),
#       scort/grampian_ifn_score.csv (冻结 IFN 分数 + Ayers6 + NK/Cyto 备用)
# 疗效标签不在导出内 → 本脚本只做表达侧; validation 脚本等标签接入
# ==============================

suppressMessages({
  library(affy); library(oligo); library(Biobase)
})
setwd("D:/GEO2R")
source("cohort_utils.R")   # symbolize(), load_signature_sets()

# ---- 使用已 RMA 的 eset (test_cel_read.R 生成, pd.clariom.s.human, 27189x223) ----
eset <- readRDS("scort/grampian_eset.rds")
y <- exprs(eset)
cat("RMA matrix:", dim(y)[1], "x", dim(y)[2], "
")

# ---- 探针 → symbol (用成功读入的平台包注释) ----
pd_pkg <- "pd.clariom.s.human"
cat("使用注释包:", pd_pkg, "\n")
con <- DBI::dbConnect(RSQLite::SQLite(),
                      system.file("extdata", paste0(pd_pkg, ".sqlite"),
                                  package = pd_pkg))
tabs <- DBI::dbListTables(con)
d <- DBI::dbGetQuery(con, "select fsetid, gene_assignment from featureSet")
ga <- d$gene_assignment
# gene_assignment 格式: "accession // symbol // desc // location"
fvec <- sapply(strsplit(ga, "//"),
               function(x) if (length(x) >= 2) trimws(x[2]) else NA_character_)
names(fvec) <- d$fsetid
DBI::dbDisconnect(con)
sym_full <- unname(fvec[colnames(y)])
cat("注释覆盖:", mean(!is.na(sym_full) & sym_full != ""), "\n")

y_sym <- symbolize(y, sym_full)
cat("symbol 化后:", dim(y_sym), "\n")
write.csv(y_sym, "scort/grampian_expr_symbol.csv")

# ---- 冻结 IFN 打分 (IFN_model_card.md 定义) ----
sig <- load_signature_sets()
sc_ifn <- GSVA::gsva(GSVA::ssgseaParam(y_sym, sig$ifn_sets), verbose = FALSE)
ifn_score <- colMeans(sc_ifn)
# Ayers 6 备用
ayers <- data.frame(HUGO.symbols = c("IFNG","STAT1","IDO1","CXCL9","CXCL10","HLA-DRA"))
# MCP-counter 官方签名 NK/Cyto 备用
genes_sig <- read.table("mcpcounter_genes_official.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
mcp_sc <- MCPcounter::MCPcounter.estimate(y_sym, featuresType = "HUGO_symbols", genes = genes_sig)

out <- data.frame(
  sample_id = names(ifn_score),
  IFN = as.numeric(ifn_score),
  IFNG_Ayers6 = as.numeric(colMeans(y_sym[intersect(c("IFNG","STAT1","IDO1","CXCL9","CXCL10","HLA-DRA"), rownames(y_sym)), ])),
  NK_mcp = as.numeric(mcp_sc["NK cells", names(ifn_score)]),
  Cyto_mcp = as.numeric(mcp_sc["Cytotoxic lymphocytes", names(ifn_score)]))
write.csv(out, "scort/grampian_ifn_score.csv", row.names = FALSE)
cat("\n===== Grampian 表达侧完成 =====\n")
cat("IFN 分数分布: min", round(min(out$IFN),3), "median", round(median(out$IFN),3),
    "max", round(max(out$IFN),3), "\n")
