# ==============================
# GSE87211 生存分析: IFN 评分 vs DFS
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(survival)
  library(pROC)
})

eset <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                  AnnotGPL = FALSE, getGPL = FALSE)$eset
y <- exprs(eset)
clin <- pData(eset)

# 临床字段
tissue_col <- grep("^tissue", colnames(clin), value = TRUE)[1]
is_tum <- clin[[tissue_col]] == "rectal tumor"
# 临床随访列含空字符串(""), as.numeric 会发 "NAs introduced by coercion";
# 空值随后由 !is.na(...) & dfs_time>0 剔除, 故显式抑制该已知良性提示。
num_clin <- function(x) suppressWarnings(as.numeric(as.character(x)))
dfs_time <- num_clin(clin[["disease free time (month):ch1"]])
dfs_event <- num_clin(clin[["cancer recurrance after surgery:ch1"]])
os_time <- num_clin(clin[["survival time (month):ch1"]])
os_event <- num_clin(clin[["death due to tumor:ch1"]])

ok <- is_tum & !is.na(dfs_time) & !is.na(dfs_event) & dfs_time > 0
cat("DFS 可用的肿瘤样本:", sum(ok), " 复发事件:", sum(dfs_event[ok]), "\n")

y_tum <- y[, ok]
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym_full <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y_tum)]
sym_first <- sapply(strsplit(sym_full, " */// *|, *"), `[`, 1)

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）
y_sym <- symbolize(y_tum, sym_first)

gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]

sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
ifn <- colMeans(sc_ifn)[colnames(y_tum)]

mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
sc_mcp <- gsva(ssgseaParam(y_sym, mcp_sets), verbose = FALSE)
nk <- sc_mcp["NK cells", colnames(y_tum)]
cyto <- sc_mcp["Cytotoxic lymphocytes", colnames(y_tum)]

srv <- data.frame(
  time = dfs_time[ok], event = dfs_event[ok],
  ifn = as.numeric(ifn), nk = as.numeric(nk), cyto = as.numeric(cyto)
)

# ---- Cox 连续变量 ----
cat("\n===== Cox 回归(连续 z 分数, DFS) =====\n")
for (v in c("ifn", "nk", "cyto")) {
  srv$zv <- as.numeric(scale(srv[[v]]))
  fit <- coxph(Surv(time, event) ~ zv, data = srv)
  s <- summary(fit)
  cat(sprintf("%-6s HR = %.3f (95%%CI %.3f-%.3f)  P = %.4f\n",
              v, s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5]))
}

# ---- KM: IFN 高 vs 低(中位数切分) ----
srv$ifn_grp <- ifelse(srv$ifn > median(srv$ifn), "IFN-high", "IFN-low")
sdiff <- survdiff(Surv(time, event) ~ ifn_grp, data = srv)
p_km <- 1 - pchisq(sdiff$chisq, 1)
cat("\n===== KM(log-rank): IFN-high vs IFN-low =====\n")
cat("log-rank P =", signif(p_km, 3), "\n")
print(table(srv$ifn_grp, srv$event))

png("GSE87211_IFN_DFS_KM.png", width = 1000, height = 850, res = 130)
km <- survfit(Surv(time, event) ~ ifn_grp, data = srv)
plot(km, col = c("tomato", "steelblue"), lwd = 2.5,
     xlab = "Months", ylab = "Disease-free survival",
     main = sprintf("GSE87211 DFS by IFN score (log-rank P = %.3f)", p_km),
     mark.time = TRUE)
legend("bottomleft", legend = c("IFN-high", "IFN-low"),
       col = c("tomato", "steelblue"), lwd = 2.5)
dev.off()

# ---- OS 顺带验证 ----
ok_os <- is_tum & !is.na(os_time) & !is.na(os_event) & os_time > 0
cat("\n===== OS 可用样本:", sum(ok_os), "=====\n")
srv_os <- data.frame(time = os_time[ok_os], event = os_event[ok_os])
srv_os$ifn <- as.numeric(ifn[colnames(y)[ok_os]])
srv_os$zv <- as.numeric(scale(srv_os$ifn))
fit_os <- coxph(Surv(time, event) ~ zv, data = srv_os)
s <- summary(fit_os)
cat(sprintf("IFN(连续) OS: HR = %.3f (%.3f-%.3f) P = %.4f\n",
            s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5]))

cat("\n===== 完成 =====\n")
