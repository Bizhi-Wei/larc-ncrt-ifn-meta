# =====================================================================
# cohort_utils.R — 跨队列共享的表达矩阵注释与特征构建工具
#
# 本文件只收口原本在各队列脚本中逐字重复的纯数据处理逻辑，不改变任何
# 统计口径或数值结果。队列特异的平台注释、疗效分组逻辑仍留在各脚本内。
# 采用与 meta_utils.R 相同的“谁用谁 source”的自包含模式；GSVA/AnnotationDbi
# 使用显式命名空间，不依赖 attach 顺序。
# =====================================================================

# ---- 1. 特征基因集 ----------------------------------------------------
# IFN = Hallmark INTERFERON_ALPHA/GAMMA_RESPONSE 两个集合（后续取均值）；
# MCP-counter 公开标志基因集，用于 ssGSEA（非 MCPcounter 反卷积算法）。
load_signature_sets <- function(gmt_file = "h.all.symbols.gmt",
                                mcp_file = "mcpcounter_genes.txt") {
  gmt <- readLines(gmt_file)
  hallmark <- lapply(strsplit(gmt, "\t"), function(x) setNames(list(x[-c(1, 2)]), x[1]))
  hallmark <- unlist(hallmark, recursive = FALSE)
  hallmark <- lapply(hallmark, function(g) g[g != ""])
  ifn_sets <- hallmark[c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                         "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
  mcp <- read.delim(mcp_file, stringsAsFactors = FALSE)
  mcp_sets <- split(mcp$HUGO.symbols, mcp$Cell.population)
  list(ifn_sets = ifn_sets, mcp_sets = mcp_sets)
}

# ---- 2. 探针/基因 → HGNC symbol --------------------------------------
# 多别名标记（" /// " 或 ", "）取第一个符号；同一 symbol 命中多个探针时，
# 按符号升序、行方差降序排序后保留方差最大的那一条。该实现是各脚本原有
# “简单版 symbolize”的超集：对本就不含分隔符的输入，别名拆分步为恒等操作。
symbolize <- function(y, symbols) {
  sym_first <- sapply(strsplit(ifelse(is.na(symbols), "", symbols), " */// *|, *"), `[`, 1)
  keep <- !is.na(sym_first) & sym_first != ""
  y <- y[keep, , drop = FALSE]; sym_first <- sym_first[keep]
  v <- apply(y, 1, var); o <- order(sym_first, -v)
  y <- y[o, , drop = FALSE]; sym_first <- sym_first[o]
  dup <- duplicated(sym_first)
  y <- y[!dup, , drop = FALSE]; rownames(y) <- sym_first[!dup]
  y
}

# ---- 3. GenBank ACCNUM → HGNC symbol（平台注释回补） -----------------
gb2sym <- function(gb) {
  gb <- sub("\\..*", "", gb)
  m <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(gb[!is.na(gb)]),
    columns = "SYMBOL", keytype = "ACCNUM"))
  m <- m[!duplicated(m$ACCNUM), ]
  unname(setNames(m$SYMBOL, m$ACCNUM)[gb])
}

# ---- 4. 仅用于运行日志的“宽容”Hedges g -------------------------------
# 允许特征缺失 / 某组不足时得到 NA 而不中断（IGF2/L1CAM 可能平台缺注释）。
# 正式 meta 效应量一律使用 meta_utils.R 的 hedges_g()，二者点估计公式相同。
.hedges_g_soft <- function(x, g) {
  x1 <- x[g == "good"]; x0 <- x[g == "poor"]
  n1 <- length(x1); n0 <- length(x0)
  sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / (n1 + n0 - 2))
  (1 - 3 / (4 * (n1 + n0) - 9)) * (mean(x1) - mean(x0)) / sp
}

# ---- 5. 由 symbol 化表达矩阵构建队列特征表 ----------------------------
# 列：IFN（两个 Hallmark IFN 集合 ssGSEA 的均值）、IGF2、L1CAM、
# NK / Cyto（MCP-counter 标志 ssGSEA）、response；剔除无疗效标签样本，
# 写出 features_<cohort>.csv 并打印一行队列摘要。
build_cohort_features <- function(y_sym, grp, cohort, ifn_sets, mcp_sets,
                                  write = TRUE, verbose = TRUE) {
  sc_ifn <- GSVA::gsva(GSVA::ssgseaParam(y_sym, ifn_sets), verbose = FALSE)
  sc_mcp <- GSVA::gsva(GSVA::ssgseaParam(y_sym, mcp_sets), verbose = FALSE)
  nk_row <- grep("NK", rownames(sc_mcp), value = TRUE)[1]
  cy_row <- grep("Cytotoxic", rownames(sc_mcp), value = TRUE)[1]
  get1 <- function(g) if (g %in% rownames(y_sym)) as.numeric(y_sym[g, ]) else rep(NA_real_, ncol(y_sym))
  feat <- data.frame(IFN = colMeans(sc_ifn),
                     IGF2 = get1("IGF2"), L1CAM = get1("L1CAM"),
                     NK = as.numeric(sc_mcp[nk_row, ]),
                     Cyto = as.numeric(sc_mcp[cy_row, ]),
                     response = grp[colnames(y_sym)],
                     row.names = colnames(y_sym))
  feat <- feat[!is.na(feat$response), ]
  if (write) write.csv(feat, paste0("features_", cohort, ".csv"))
  if (verbose) {
    cat(sprintf("[%s] n=%d (good %d / poor %d)  g: IFN=%.2f NK=%.2f Cyto=%.2f IGF2=%.2f L1CAM=%.2f\n",
                cohort, nrow(feat), sum(feat$response == "good"), sum(feat$response == "poor"),
                .hedges_g_soft(feat$IFN, feat$response), .hedges_g_soft(feat$NK, feat$response),
                .hedges_g_soft(feat$Cyto, feat$response),
                .hedges_g_soft(feat$IGF2, feat$response), .hedges_g_soft(feat$L1CAM, feat$response)))
  }
  invisible(feat)
}
