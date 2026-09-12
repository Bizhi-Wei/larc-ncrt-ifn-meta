# ==============================
# 四个初始队列 Hallmark 通路 camera 富集分析
# 方向统一: Up = 疗效好组(CR/pCR/Responder)高表达
# 注：脚本和旧输出文件名保留 "3cohorts" 仅为向后兼容。
# ==============================

suppressMessages({
  library(GEOquery)
  library(limma)
  library(DESeq2)
  library(ggplot2)
})

# ---- 读取 Hallmark 基因集 ----
gmt <- readLines("h.all.symbols.gmt")
hallmark <- lapply(strsplit(gmt, "\t"), function(x) {
  setNames(list(x[-c(1, 2)]), x[1])
})
hallmark <- unlist(hallmark, recursive = FALSE)
hallmark <- lapply(hallmark, function(g) g[g != ""])
cat("Hallmark 基因集数:", length(hallmark), "\n")

# 通用 camera 函数: y=表达矩阵(行=基因符号), group 两水平(好/差)
run_camera <- function(y, symbols, group, good_level) {
  group <- factor(group)
  group <- relevel(group, ref = setdiff(levels(group), good_level))
  design <- model.matrix(~ group)
  # 同一符号多行时取方差最大的一行
  keep <- !is.na(symbols) & symbols != ""
  y <- y[keep, ]; symbols <- symbols[keep]
  v <- apply(y, 1, var)
  y <- y[order(symbols, -v), ]; symbols <- symbols[order(symbols, -v)]
  dup <- duplicated(symbols)
  y <- y[!dup, ]; symbols <- symbols[!dup]
  rownames(y) <- symbols
  idx <- ids2indices(hallmark, symbols)
  cam <- camera(y, idx, design, contrast = 2)
  cam$set <- rownames(cam)
  cam[, c("set", "NGenes", "Direction", "PValue", "FDR")]
}

# ---- 队列1: GSE209746 (VST, CR vs iCR) ----
vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
y1 <- assay(vsd)
s8a <- read.csv("paper_S8A_results.csv", stringsAsFactors = FALSE)
ens2sym <- setNames(s8a$hgnc_symbol, s8a$ensembl_gene_id)
sym1 <- ens2sym[sub("\\..*$", "", rownames(y1))]
grp1 <- labels$response[match(colnames(y1), labels$patient)]
grp1 <- ifelse(grp1 == "CR", "good", "poor")

r1 <- run_camera(y1, sym1, grp1, "good")

# ---- 队列2: GSE35452 (limma输入矩阵, Responder vs Non-responder) ----
eset2 <- readRDS("GSE35452_eset.rds")
y2 <- exprs(eset2)
fd2 <- fData(eset2)
clin2 <- pData(eset2)
rt <- as.character(clin2$characteristics_ch1.1)
grp2 <- ifelse(grepl("Non-responder", rt, ignore.case = TRUE), "poor", "good")
sym2 <- sapply(strsplit(as.character(fd2[["Gene symbol"]]), " */// *"), `[`, 1)

r2 <- run_camera(y2, sym2, grp2, "good")

# ---- 队列3: GSE87211 (pCR vs non-pCR, 仅肿瘤样本) ----
eset3 <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
y3 <- exprs(eset3)
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
y3 <- y3[, ok3]

pmap <- read.delim("GPL13497_probe_map.tsv", stringsAsFactors = FALSE)
sym3 <- setNames(pmap$GENE_SYMBOL, pmap$ID)[rownames(y3)]

r3 <- run_camera(y3, sym3, grp3, "good")

# ---- 队列4: GSE150082 (Good vs Poor) ----
# 该表由流水线上一步 GSE150082_validation.R 用相同 camera 口径生成。
cam4 <- read.csv("GSE150082_hallmark_camera.csv", stringsAsFactors = FALSE)
r4 <- cam4[, c("set", "NGenes", "Direction", "PValue", "FDR")]

# ---- 合并三队列 ----
colnames(r1) <- c("set", "NGenes_GSE209746", "Direction_GSE209746",
                  "PValue_GSE209746", "FDR_GSE209746")
colnames(r2) <- c("set", "NGenes_GSE35452", "Direction_GSE35452",
                  "PValue_GSE35452", "FDR_GSE35452")
colnames(r3) <- c("set", "NGenes_GSE87211", "Direction_GSE87211",
                  "PValue_GSE87211", "FDR_GSE87211")
colnames(r4) <- c("set", "NGenes_GSE150082", "Direction_GSE150082",
                  "PValue_GSE150082", "FDR_GSE150082")
m <- merge(r1, r2, by = "set")
m <- merge(m, r3, by = "set")
m <- merge(m, r4, by = "set")

# 方向打分: Up=+1, Down=-1
dir_score <- function(d) ifelse(d == "Up", 1, -1)
m$score1 <- dir_score(m$Direction_GSE209746)
m$score2 <- dir_score(m$Direction_GSE35452)
m$score3 <- dir_score(m$Direction_GSE87211)
m$score4 <- dir_score(m$Direction_GSE150082)
m$consistent <- (m$score1 == m$score2) & (m$score2 == m$score3) &
  (m$score3 == m$score4)
m$consistency_score <- m$score1 + m$score2 + m$score3 + m$score4
m$min_P <- pmin(m$PValue_GSE209746, m$PValue_GSE35452,
                m$PValue_GSE87211, m$PValue_GSE150082)
m$prodP_rank <- rank(m$PValue_GSE209746 * m$PValue_GSE35452 *
                     m$PValue_GSE87211 * m$PValue_GSE150082)

# 按一致性+P值排序
m <- m[order(!m$consistent, m$PValue_GSE209746 *
             m$PValue_GSE35452 * m$PValue_GSE87211 *
             m$PValue_GSE150082), ]

show_cols <- c("set", "Direction_GSE209746", "PValue_GSE209746",
               "Direction_GSE35452", "PValue_GSE35452",
               "Direction_GSE87211", "PValue_GSE87211",
               "Direction_GSE150082", "PValue_GSE150082", "consistent")
cat("\n===== Hallmark 通路四队列结果(按一致性+P乘积排序) =====\n")
out <- m[, show_cols]
out$PValue_GSE209746 <- signif(out$PValue_GSE209746, 3)
out$PValue_GSE35452 <- signif(out$PValue_GSE35452, 3)
out$PValue_GSE87211 <- signif(out$PValue_GSE87211, 3)
out$PValue_GSE150082 <- signif(out$PValue_GSE150082, 3)
print(out, row.names = FALSE)

write.csv(m, "hallmark_3cohorts_camera.csv", row.names = FALSE)
write.csv(m, "hallmark_4cohorts_camera.csv", row.names = FALSE)

# ---- Figure S1: selected cross-cohort Hallmark pathways ----
selected_sets <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA"
)
selected_sets <- selected_sets[selected_sets %in% m$set]
plot_rows <- lapply(c("GSE209746", "GSE35452", "GSE87211", "GSE150082"), function(cohort) {
  direction <- m[[paste0("Direction_", cohort)]]
  p_value <- m[[paste0("PValue_", cohort)]]
  data.frame(
    set = m$set,
    cohort = cohort,
    signed_log10_p = ifelse(direction == "Up", 1, -1) * -log10(pmax(p_value, 1e-300)),
    label = sprintf("%s\nP=%s", ifelse(direction == "Up", "Up", "Down"),
                    format.pval(p_value, digits = 2, eps = 1e-99)),
    stringsAsFactors = FALSE
  )
})
plot_df <- do.call(rbind, plot_rows)
plot_df <- plot_df[plot_df$set %in% selected_sets, ]
pretty_set <- gsub("^HALLMARK_", "", selected_sets)
pretty_set <- gsub("_", " ", pretty_set)
plot_df$set <- factor(plot_df$set, levels = rev(selected_sets), labels = rev(pretty_set))
plot_df$cohort <- factor(plot_df$cohort,
                         levels = c("GSE209746", "GSE35452", "GSE87211", "GSE150082"))
fill_limit <- max(abs(plot_df$signed_log10_p), na.rm = TRUE)

p_s1 <- ggplot(plot_df, aes(x = cohort, y = set, fill = signed_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.45, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#D55E00", mid = "white", high = "#0072B2", midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    name = "Signed -log10(P)\nUp in responders > 0"
  ) +
  labs(x = NULL, y = NULL,
       title = "Hallmark pathway consistency across four initial cohorts") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black", size = 7.6),
    panel.grid = element_blank(),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  )

ggsave("hallmark_3cohorts_barplot.png", p_s1, width = 183, height = 150,
       units = "mm", dpi = 300, bg = "white")
ggsave("hallmark_3cohorts_barplot.pdf", p_s1, width = 183, height = 150,
       units = "mm", device = cairo_pdf, family = "Arial")
svglite::svglite("hallmark_3cohorts_barplot.svg", width = 183 / 25.4,
                 height = 150 / 25.4)
print(p_s1)
dev.off()

# ---- 关注主题速览 ----
cat("\n===== 关注主题速览 =====\n")
themes <- grep("HYPOXIA|DNA_REPAIR|E2F|MYC|OXIDATIVE|GLYCOLYSIS|MTORC1|INTERFERON|INFLAMMATORY|COMPLEMENT|IL6|IL2|TNFA|EPITHELIAL|ANGIOGENESIS|APICAL|TGFB|FATTY|CHOLESTEROL|KRAS|P53|G2M|ALLOGRAFT",
               m$set, value = TRUE)
print(m[m$set %in% themes, show_cols], row.names = FALSE)

cat("\n===== 完成 =====\n")
