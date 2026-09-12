# ==============================
# 肿瘤纯度/间质校正: IFN-疗效关联是否独立于间质含量
# 校正因子: MCP-counter Fibroblasts + Endothelial cells (ESTIMATE 包未装, 用MCP间质代理)
# 四队列: GSE209746 / GSE35452 / GSE87211 / GSE150082
# ==============================

suppressMessages({
  library(GEOquery)
  library(GSVA)
  library(DESeq2)
})

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）

gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
stroma_sets <- mcp_sets[c("Fibroblasts", "Endothelial cells")]

build <- function(y_sym, grp) {
  sc_ifn <- gsva(ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
  sc_str <- gsva(ssgseaParam(y_sym, stroma_sets), verbose = FALSE)
  data.frame(
    ifn   = as.numeric(colMeans(sc_ifn)),
    fibro = as.numeric(sc_str["Fibroblasts", ]),
    endo  = as.numeric(sc_str["Endothelial cells", ]),
    response = grp[colnames(y_sym)],
    row.names = colnames(y_sym)
  )
}

# ---- 队列1: GSE209746 ----
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym1 <- ens2sym[sub("\\..*$", "", rownames(assay(vsd)))]
y1 <- symbolize(assay(vsd), sym1)
grp1 <- ifelse(labels$response[match(colnames(y1), labels$patient)] == "CR", "good", "poor")
names(grp1) <- colnames(y1)
f1 <- build(y1, grp1)

# ---- 队列2: GSE35452 ----
eset2 <- readRDS("GSE35452_eset.rds")
fd2 <- fData(eset2)
rt <- as.character(pData(eset2)$characteristics_ch1.1)
grp2 <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
names(grp2) <- colnames(exprs(eset2))
sym2 <- sapply(strsplit(as.character(fd2[["Gene symbol"]]), " */// *"), `[`, 1)
f2 <- build(symbolize(exprs(eset2), sym2), grp2)

# ---- 队列3: GSE87211 (仅肿瘤) ----
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
names(grp3) <- colnames(y3raw)[ok3]
pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym3 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y3raw)]
sym3 <- sapply(strsplit(sym3, " */// *|, *"), `[`, 1)
f3 <- build(symbolize(y3raw[, ok3], sym3), grp3)

# ---- 队列4: GSE150082 ----
eset4 <- GEOquery:::parseGSEMatrix("GSE150082_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
y4 <- exprs(eset4)
grp4 <- ifelse(pData(eset4)$`response:ch1` == "Good", "good", "poor")
names(grp4) <- colnames(y4)
sym4 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y4)]
sym4 <- sapply(strsplit(sym4, " */// *|, *"), `[`, 1)
f4 <- build(symbolize(y4, sym4), grp4)

# ==============================
# 每队列: 间质 vs 疗效 / IFN~response+间质 / 分层
# ==============================
detail_rows <- list()

analyze <- function(f, tag) {
  f$ifn_z <- as.numeric(scale(f$ifn))
  f$fibro_z <- as.numeric(scale(f$fibro))
  f$endo_z <- as.numeric(scale(f$endo))
  f$y <- ifelse(f$response == "good", 1, 0)

  cat(sprintf("\n===== %s (n=%d, good=%d) =====\n", tag, nrow(f), sum(f$y)))
  # 间质本身与疗效的关系
  for (v in c("fibro", "endo")) {
    tt <- t.test(f[[v]][f$response == "good"], f[[v]][f$response == "poor"])
    cat(sprintf("  %-6s 好组均值=%.3f 差组均值=%.3f  P=%.4f (%s)\n", v,
                mean(f[[v]][f$response == "good"]),
                mean(f[[v]][f$response == "poor"]), tt$p.value,
                ifelse(mean(f[[v]][f$response == "good"]) >
                       mean(f[[v]][f$response == "poor"]), "good高", "poor高")))
  }
  # 单变量 vs 校正间质后: 线性模型 IFN ~ response (+fibro+endo)
  m0 <- lm(ifn_z ~ y, data = f)
  m1 <- lm(ifn_z ~ y + fibro_z + endo_z, data = f)
  cat(sprintf("  IFN~response:          beta=%.3f P=%.4f\n",
              coef(summary(m0))["y", 1], coef(summary(m0))["y", 4]))
  cat(sprintf("  IFN~response+间质校正: beta=%.3f P=%.4f\n",
              coef(summary(m1))["y", 1], coef(summary(m1))["y", 4]))
  # 反向: logistic response ~ IFN (+间质)
  g0 <- glm(y ~ ifn_z, data = f, family = binomial)
  g1 <- glm(y ~ ifn_z + fibro_z + endo_z, data = f, family = binomial)
  cat(sprintf("  response~IFN:          OR=%.2f P=%.4f\n",
              exp(coef(g0)["ifn_z"]), coef(summary(g0))["ifn_z", 4]))
  cat(sprintf("  response~IFN+间质校正: OR=%.2f P=%.4f\n",
              exp(coef(g1)["ifn_z"]), coef(summary(g1))["ifn_z", 4]))
  # IFN 与间质的相关性(信号是否只是纯度反映)
  cat(sprintf("  cor(IFN, Fibro) = %.3f   cor(IFN, Endo) = %.3f\n",
              cor(f$ifn, f$fibro), cor(f$ifn, f$endo)))

  linear_row <- function(fit, model_name, adjustment) {
    ci <- confint.default(fit)["y", ]
    data.frame(
      section = "Stromal adjustment", cohort = tag,
      population = "Pretreatment tumors with response labels",
      outcome = "IFN score (z)", model = model_name,
      focal_term = "Response (good vs poor)", effect_measure = "Beta",
      estimate = unname(coef(fit)["y"]), ci_lo = unname(ci[1]),
      ci_hi = unname(ci[2]), p_value = coef(summary(fit))["y", 4],
      n = nobs(fit), events = sum(model.frame(fit)$y == 1), adjustment = adjustment,
      note = "Positive beta indicates higher IFN score in responders.",
      source = "tumor_purity_adjust_models.csv",
      stringsAsFactors = FALSE
    )
  }
  logistic_row <- function(fit, model_name, adjustment) {
    ci <- exp(confint.default(fit)["ifn_z", ])
    data.frame(
      section = "Stromal adjustment", cohort = tag,
      population = "Pretreatment tumors with response labels",
      outcome = "Good response", model = model_name,
      focal_term = "IFN score (per SD)", effect_measure = "Odds ratio",
      estimate = unname(exp(coef(fit)["ifn_z"])), ci_lo = unname(ci[1]),
      ci_hi = unname(ci[2]), p_value = coef(summary(fit))["ifn_z", 4],
      n = nobs(fit), events = sum(model.frame(fit)$y == 1), adjustment = adjustment,
      note = "OR above 1 indicates greater odds of good response.",
      source = "tumor_purity_adjust_models.csv",
      stringsAsFactors = FALSE
    )
  }
  detail_rows[[tag]] <<- rbind(
    linear_row(m0, "Linear regression", "None"),
    linear_row(m1, "Linear regression", "Fibroblast and endothelial ssGSEA scores"),
    logistic_row(g0, "Logistic regression", "None"),
    logistic_row(g1, "Logistic regression", "Fibroblast and endothelial ssGSEA scores")
  )

  data.frame(cohort = tag,
             beta_unadj = coef(summary(m0))["y", 1], P_unadj = coef(summary(m0))["y", 4],
             beta_adj = coef(summary(m1))["y", 1], P_adj = coef(summary(m1))["y", 4],
             OR_unadj = exp(coef(g0)["ifn_z"]), OR_adj = exp(coef(g1)["ifn_z"]),
             P_logit_adj = coef(summary(g1))["ifn_z", 4],
             cor_ifn_fibro = cor(f$ifn, f$fibro), cor_ifn_endo = cor(f$ifn, f$endo))
}

res <- rbind(analyze(f1, "GSE209746"), analyze(f2, "GSE35452"),
             analyze(f3, "GSE87211"), analyze(f4, "GSE150082"))
res$beta_unadj <- round(res$beta_unadj, 3); res$beta_adj <- round(res$beta_adj, 3)
res$P_unadj <- signif(res$P_unadj, 3); res$P_adj <- signif(res$P_adj, 3)
res$OR_unadj <- round(res$OR_unadj, 2); res$OR_adj <- round(res$OR_adj, 2)
res$P_logit_adj <- signif(res$P_logit_adj, 3)
res$cor_ifn_fibro <- round(res$cor_ifn_fibro, 3)
res$cor_ifn_endo <- round(res$cor_ifn_endo, 3)
cat("\n===== 汇总 =====\n")
print(res, row.names = FALSE)
write.csv(res, "tumor_purity_adjust.csv", row.names = FALSE)
write.csv(do.call(rbind, detail_rows), "tumor_purity_adjust_models.csv", row.names = FALSE)

cat("\n===== 完成 =====\n")
