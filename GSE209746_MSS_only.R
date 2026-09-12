# ==============================
# GSE209746 MSS-only 敏感性分析
# 剔除 MSI/NA 样本后重验: IFN签名 / IGF2/L1CAM / 免疫浸润 / 免疫热聚类 / 15基因DESeq2
# ==============================

suppressMessages({
  library(DESeq2)
  library(GSVA)
  library(pROC)
})

# ---- 数据加载 ----
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
s1d <- read.csv("paper_S1D_bindea_scores.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)

source("cohort_utils.R")  # symbolize() 等跨队列公共工具（原本地定义已收口）

y1 <- symbolize(assay(vsd), ens2sym[sub("\\..*$", "", rownames(assay(vsd)))])
grp1 <- ifelse(labels$response[match(colnames(y1), labels$patient)] == "CR", "good", "poor")
names(grp1) <- colnames(y1)

msi_raw <- setNames(s1d$MSI_status, s1d$patient)[colnames(y1)]
cat("MSI_status 原始取值:", paste(names(table(msi_raw, useNA = "ifany")),
                                  table(msi_raw, useNA = "ifany"), collapse = " / "), "\n")
# 0/1 编码: 0=MSS, 1=MSI (与 GEO series matrix 字符串标注计数一致: MSI=5)
msi <- ifelse(msi_raw == 1, "MSI", ifelse(msi_raw == 0, "MSS", NA))
names(msi) <- colnames(y1)

# ---- MSI × response 交叉表(全部93例) ----
cat("\n===== MSI × 疗效(全部", ncol(y1), "例) =====\n")
print(table(msi = msi, response = grp1, useNA = "ifany"))
f <- fisher.test(table(msi = msi, response = grp1, useNA = "no"))
cat(sprintf("MSI 与疗效关联 Fisher P = %.4f (OR = %.2f)\n", f$p.value, f$estimate))

# ---- MSS-only 子集 ----
mss <- colnames(y1)[msi == "MSS" & !is.na(msi)]
cat("\nMSS-only 样本数:", length(mss), " CR:", sum(grp1[mss] == "good"),
    " iCR:", sum(grp1[mss] == "poor"), "\n")
grp_mss <- grp1[mss]

# ---- 基因集 ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
mcp <- read.delim("mcpcounter_genes.txt", stringsAsFactors = FALSE)
mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)

hedges_g <- function(x1, x2) {  # x1=good, x2=poor
  n1 <- length(x1); n2 <- length(x2)
  sp <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
  d <- (mean(x1) - mean(x2)) / sp
  (1 - 3 / (4 * (n1 + n2) - 9)) * d
}

# ==============================
# A. 全队列 vs MSS-only: IFN / IGF2 / L1CAM / NK / Cyto 对照
# ==============================
run_block <- function(sams, tag) {
  yy <- y1[, sams]; gg <- grp1[sams]
  sc_ifn <- gsva(ssgseaParam(yy, ifn_sets), verbose = FALSE)
  ifn <- colMeans(sc_ifn)
  sc_mcp <- gsva(ssgseaParam(yy, mcp_sets), verbose = FALSE)
  feats <- list(
    IFN   = ifn,
    IGF2  = yy["IGF2", ],
    L1CAM = yy["L1CAM", ],
    NK    = sc_mcp["NK cells", ],
    Cyto  = sc_mcp["Cytotoxic lymphocytes", ]
  )
  out <- data.frame()
  for (nm in names(feats)) {
    x <- feats[[nm]]
    x1 <- x[gg == "good"]; x2 <- x[gg == "poor"]
    tt <- t.test(x1, x2)
    r <- roc(gg, as.numeric(x), levels = c("poor", "good"), direction = "<", quiet = TRUE)
    out <- rbind(out, data.frame(
      feature = nm, cohort = tag,
      good_mean = round(mean(x1), 3), poor_mean = round(mean(x2), 3),
      hedges_g = round(hedges_g(x1, x2), 3),
      P = signif(tt$p.value, 3), AUC = round(as.numeric(auc(r)), 3)
    ))
  }
  list(res = out, sc_mcp = sc_mcp)
}

rb_full <- run_block(colnames(y1), "full_n93")
rb_mss  <- run_block(mss, "MSS_only")
cmp <- rbind(rb_full$res, rb_mss$res)
cat("\n===== 全队列 vs MSS-only 特征对照 =====\n")
print(cmp, row.names = FALSE)
write.csv(cmp, "GSE209746_MSS_only_features.csv", row.names = FALSE)

g_ci <- function(g, n1, n2) {
  se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2 - 2)))
  c(g - qnorm(0.975) * se, g + qnorm(0.975) * se)
}
ifn_full <- cmp[cmp$feature == "IFN" & cmp$cohort == "full_n93", ]
ifn_mss <- cmp[cmp$feature == "IFN" & cmp$cohort == "MSS_only", ]
ci_full <- g_ci(ifn_full$hedges_g, sum(grp1 == "good"), sum(grp1 == "poor"))
ci_mss <- g_ci(ifn_mss$hedges_g, sum(grp_mss == "good"), sum(grp_mss == "poor"))
msi_detail <- rbind(
  data.frame(
    section = "MSI sensitivity", cohort = "GSE209746",
    population = "Tumors with nonmissing MSI status", outcome = "Good response",
    model = "Fisher exact test", focal_term = "MSI vs MSS", effect_measure = "Odds ratio",
    estimate = unname(f$estimate), ci_lo = unname(f$conf.int[1]), ci_hi = unname(f$conf.int[2]),
    p_value = f$p.value, n = sum(!is.na(msi)), events = sum(grp1[!is.na(msi)] == "good"),
    adjustment = "None",
    note = "Three tumors with missing MSI status were excluded.",
    source = "GSE209746_MSI_sensitivity.csv", stringsAsFactors = FALSE
  ),
  data.frame(
    section = "MSI sensitivity", cohort = "GSE209746",
    population = "All response-labelled tumors", outcome = "IFN score",
    model = "Standardized mean difference", focal_term = "Good vs poor response",
    effect_measure = "Hedges g", estimate = ifn_full$hedges_g,
    ci_lo = ci_full[1], ci_hi = ci_full[2], p_value = ifn_full$P,
    n = ncol(y1), events = sum(grp1 == "good"), adjustment = "None",
    note = "Positive g indicates higher IFN score in responders.",
    source = "GSE209746_MSI_sensitivity.csv", stringsAsFactors = FALSE
  ),
  data.frame(
    section = "MSI sensitivity", cohort = "GSE209746",
    population = "MSS-only tumors", outcome = "IFN score",
    model = "Standardized mean difference", focal_term = "Good vs poor response",
    effect_measure = "Hedges g", estimate = ifn_mss$hedges_g,
    ci_lo = ci_mss[1], ci_hi = ci_mss[2], p_value = ifn_mss$P,
    n = length(mss), events = sum(grp_mss == "good"), adjustment = "MSS-only restriction",
    note = "MSI and MSI-missing tumors were excluded.",
    source = "GSE209746_MSI_sensitivity.csv", stringsAsFactors = FALSE
  )
)
write.csv(msi_detail, "GSE209746_MSI_sensitivity.csv", row.names = FALSE)

# ==============================
# B. 免疫热聚类复现 + MSS-only 重算
# ==============================
cluster_immune <- function(sc_mcp, grp, tag) {
  z <- t(scale(t(sc_mcp)))                    # 按细胞类型 z 化
  hc <- hclust(dist(t(z)), method = "ward.D2")
  cl <- cutree(hc, k = 4)
  ord <- order(sapply(split(colMeans(z), cl), mean))  # 冷→热排序
  ig <- paste0("IG", ord[cl])
  names(ig) <- colnames(z)
  cr_hot <- mean(grp[names(ig)][ig %in% c("IG3", "IG4")] == "good")
  cr_cold <- mean(grp[names(ig)][ig %in% c("IG1", "IG2")] == "good")
  tab <- table(immune = ifelse(ig %in% c("IG3", "IG4"), "hot_IG34", "cold_IG12"),
               response = grp[names(ig)])
  ft <- fisher.test(tab)
  cat(sprintf("\n[%s] 聚类分布: %s\n", tag, paste(names(table(ig)), table(ig), collapse = " ")))
  cat(sprintf("[%s] CR率 IG3/4 = %.0f%% (%d/%d) vs IG1/2 = %.0f%% (%d/%d), Fisher P = %.4f\n",
              tag, cr_hot * 100, sum(grp[names(ig)][ig %in% c("IG3","IG4")] == "good"),
              sum(ig %in% c("IG3","IG4")),
              cr_cold * 100, sum(grp[names(ig)][ig %in% c("IG1","IG2")] == "good"),
              sum(ig %in% c("IG1","IG2")), ft$p.value))
  invisible(list(ig = ig, p = ft$p.value, tab = tab))
}

cat("\n===== 免疫热聚类: 全队列(复现验证) vs MSS-only =====\n")
ic_full <- cluster_immune(rb_full$sc_mcp, grp1[colnames(rb_full$sc_mcp)], "全队列93例")
ic_mss  <- cluster_immune(rb_mss$sc_mcp, grp_mss[colnames(rb_mss$sc_mcp)], "MSS-only")

# ==============================
# C. DESeq2 MSS-only: 15 基因
# ==============================
cat("\n===== DESeq2 MSS-only(15基因) =====\n")
dds0 <- readRDS("GSE209746_dds_loose.rds")
cts <- counts(dds0)
cd <- as.data.frame(colData(dds0))
cd$patient <- rownames(cd)
cd$resp <- labels$response[match(cd$patient, labels$patient)]
cd$msi <- msi[cd$patient]
keep_s <- !is.na(cd$msi) & cd$msi == "MSS" & !is.na(cd$resp)
# 进入设计前先把 resp 显式 factor 并锁定 CR 为参照水平,
# 避免 DESeq2 "design variables are characters, converting to factors" 提示。
cd_mss <- cd[keep_s, ]
cd_mss$resp <- relevel(factor(cd_mss$resp, levels = c("CR", "iCR")), ref = "CR")
dds <- DESeqDataSetFromMatrix(cts[, keep_s], cd_mss, ~ resp)
keep_g <- rowSums(counts(dds)) >= 10
dds <- DESeq(dds[keep_g, ], quiet = TRUE)
res <- results(dds, contrast = c("resp", "iCR", "CR"))
res$symbol <- ens2sym[sub("\\..*$", "", rownames(res))]

genes15 <- data.frame(
  symbol = c("VAX1","REG3A","IGF2","NHLH2","SOX2","KLK6","L1CAM","INHBB",
             "CPS1","ALOX12B","TDRD1","UNC13A","TENM1","ADCY2","ZNF300"),
  paper_dir = c(rep("iCR_high", 8), rep("CR_high", 7)))
tab15 <- data.frame()
for (i in seq_len(nrow(genes15))) {
  g <- genes15$symbol[i]
  hit <- which(res$symbol == g & !is.na(res$symbol))
  if (length(hit) == 0) {
    tab15 <- rbind(tab15, data.frame(symbol = g, paper_dir = genes15$paper_dir[i],
                                     log2FC = NA, padj = NA, note = "filtered"))
    next
  }
  hit <- hit[which.min(res$pvalue[hit])]
  d <- ifelse(res$log2FoldChange[hit] > 0, "iCR_high", "CR_high")
  tab15 <- rbind(tab15, data.frame(
    symbol = g, paper_dir = genes15$paper_dir[i],
    log2FC = round(res$log2FoldChange[hit], 2),
    padj = signif(res$padj[hit], 3),
    note = ifelse(d == genes15$paper_dir[i], "dir_same", "dir_opposite")))
}
print(tab15, row.names = FALSE)
sig_same <- sum(tab15$note == "dir_same" & !is.na(tab15$padj) & tab15$padj < 0.05 &
                abs(tab15$log2FC) > 1)
cat(sprintf("MSS-only 中 |log2FC|>1 且 FDR<0.05 且方向一致: %d/15\n", sig_same))
write.csv(tab15, "GSE209746_MSS_only_15genes.csv", row.names = FALSE)

cat("\n===== 完成 =====\n")
