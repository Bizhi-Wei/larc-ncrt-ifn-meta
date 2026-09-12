# ==============================
# 主图绘制 (ggplot2, 发表级)
# A: IFN signature nCRT 队列森林图
# B: NK + Cytotoxic nCRT 队列森林图 (左右两幅)
# C: GSE209746 免疫分组(IG1-4)CR率
# D: GSE87211 DFS by pCR (pCR 零复发)
# 输出: 单面板 PNG + 2x2 组合主图 Fig1_main.png (300 dpi)
# ==============================

suppressMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(GEOquery)
  library(survival)
  library(scales)
})

TXT <- "sans"
source("meta_utils.R")
fmt_p <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.4f", p))
# 避免 sprintf 把 -0.0049 印成 "-0.00"
fmt_g <- function(x) { x[abs(x) < 0.005] <- 0; sprintf("%.2f", x) }

# ---------- 森林图: 显式 y 坐标 ----------
# df: cohort, g, se, lo, hi
forest_gg <- function(df, title = NULL, show_ylab = TRUE,
                      xlab = "Hedges g (> 0: higher in responders)",
                      pool_override = NULL) {
  # Pooled 统计默认由逐队列值重算；主图传入落盘的 meta_pooled_results.csv
  # 结果(pool_override)，使图与结果表严格同源，避免绘图时重复 REML 优化。
  pool <- if (is.null(pool_override)) pool_random_effects(df$g, df$se) else pool_override
  n <- nrow(df)
  df$y <- n - seq_len(n) + 1.8          # 第1队列在最上
  pool_y <- 1                            # 合并行在底部
  hdr_y <- n + 1.55                      # 表头
  het_y <- 0.02
  x_txt <- 2.05                          # 右侧文字列
  x_min <- min(-0.75, min(df$lo) - 0.08, pool$pi_lo - 0.10)
  x_max <- max(1.85, pool$pi_hi + 0.10)

  dd_text <- data.frame(
    y = c(df$y, pool_y),
    lab = c(sprintf("%s (%s, %s)", fmt_g(df$g), fmt_g(df$lo), fmt_g(df$hi)),
            sprintf("%s (%s, %s), P=%s", fmt_g(pool$g), fmt_g(pool$lo), fmt_g(pool$hi), fmt_p(pool$p))),
    bold = c(rep(FALSE, n), TRUE)
  )
  y_br <- c(df$y, pool_y)
  y_lab <- c(df$cohort, "Pooled (REML-mKH)")
  # 仅合并行刻度标签加粗: 用 plotmath expression 向量实现(普通标签为文本、
  # 末行 bold), 避免向量化 element_text(face=) 的非官方用法。
  y_disp <- c(as.expression(lapply(y_lab[seq_len(n)], function(s) bquote(.(s)))),
              expression(bold("Pooled (REML-mKH)")))

  g <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "gray55", linewidth = 0.4) +
    geom_vline(xintercept = pool$g, linetype = 3, color = "#D55E00", alpha = 0.65,
               linewidth = 0.4) +
    geom_errorbar(data = df, aes(y = y, xmin = lo, xmax = hi),
                  orientation = "y", width = 0.16, linewidth = 0.55) +
    geom_point(data = df, aes(y = y, x = g, size = 1 / se^2),
               shape = 15, color = "black") +
    geom_polygon(data = data.frame(x = c(pool$lo, pool$g, pool$hi, pool$g),
                                   y = pool_y + c(0, 0.24, 0, -0.24)),
                 aes(x = x, y = y), fill = "#D55E00") +
    # 95% 预测区间: 菱形下方的细横杠 + 端帽
    annotate("segment", x = pool$pi_lo, xend = pool$pi_hi,
             y = pool_y - 0.38, yend = pool_y - 0.38,
             color = "#D55E00", linewidth = 0.7) +
    annotate("segment", x = pool$pi_lo, xend = pool$pi_lo,
             y = pool_y - 0.30, yend = pool_y - 0.46,
             color = "#D55E00", linewidth = 0.7) +
    annotate("segment", x = pool$pi_hi, xend = pool$pi_hi,
             y = pool_y - 0.30, yend = pool_y - 0.46,
             color = "#D55E00", linewidth = 0.7) +
    annotate("text", x = pool$pi_hi, y = pool_y - 0.38, hjust = -0.12,
             label = sprintf("95%% PI: %.2f to %.2f", pool$pi_lo, pool$pi_hi),
             size = 2.7, color = "#D55E00", family = TXT) +
    geom_text(data = dd_text, aes(x = x_txt, y = y, label = lab,
                                  fontface = ifelse(bold, "bold", "plain")),
              hjust = 0, size = 3, family = TXT) +
    annotate("text", x = x_txt, y = hdr_y, label = "g (95% CI)",
             fontface = "bold", hjust = 0, size = 3.1, family = TXT) +
    annotate("text", x = pool$g, y = het_y,
             label = sprintf("I2 = %.0f%%, Q P = %.3f", pool$I2, pool$Qp),
             size = 2.9, color = "gray30", family = TXT) +
    scale_size_continuous(range = c(1.2, 3.2), guide = "none") +
    scale_x_continuous(breaks = seq(-0.5, 1.5, 0.5)) +
    scale_y_continuous(breaks = y_br, labels = y_disp, expand = c(0, 0.15)) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.18, hdr_y + 0.45),
                    clip = "off") +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 11, base_family = TXT) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
          axis.text.y = element_text(size = 9.5, color = "black"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(6, 150, 6, 6))
  # 合并行 y 标签加粗已通过 scale_y_continuous(labels = y_disp) 的 plotmath 实现
  g <- g + theme(axis.text.y = element_text(size = 9.5, color = "black"))
  if (!show_ylab) {
    g <- g + theme(axis.text.y = element_blank())
  }
  g
}

# ---------- 读取 meta 数据 ----------
meta <- read.csv("meta_all_cohorts.csv", stringsAsFactors = FALSE)
meta <- meta[meta$treatment_scope == "nCRT", ]
# 主分析(9 nCRT)的落盘合并结果，作为森林图 Pooled 行的唯一数据源
pooled_tab <- read.csv("meta_pooled_results.csv", stringsAsFactors = FALSE)
pooled_tab <- pooled_tab[pooled_tab$analysis == "primary_ncrt9", ]
pool_from_disk <- function(ft) {
  r <- pooled_tab[pooled_tab$feature == ft, ]
  stopifnot(nrow(r) == 1L)
  list(g = r$estimate, lo = r$ci_lo, hi = r$ci_hi, p = r$p, I2 = r$I2,
       Qp = r$Qp, pi_lo = r$pi_lo, pi_hi = r$pi_hi)
}
label_map <- c(
  GSE209746 = "GSE209746 (RNA-seq, CR/iCR",
  GSE35452  = "GSE35452 (Affymetrix, TRG",
  GSE87211  = "GSE87211 (Agilent, ypT0N0",
  GSE150082 = "GSE150082 (Agilent 2col, TRG",
  GSE45404  = "GSE45404 (Affymetrix, Mandard",
  GSE119409 = "GSE119409 (Affymetrix, sens/resist",
  GSE133057 = "GSE133057 (Illumina, AJCC0-1/2-3",
  GSE53781  = "GSE53781 (CodeLink, TRG1-2/3-5",
  GSE94104  = "GSE94104 (Illumina-FFPE, TRG3/1-2"
)
meta$cohort <- paste0(label_map[meta$cohort], ", n=", meta$n, ")")
mk_df <- function(feat) meta[meta$feature == feat, c("cohort", "g", "se", "lo", "hi")]

n_primary <- length(unique(meta$cohort))
pA <- forest_gg(mk_df("IFN"), title = sprintf("A  IFN signature (%d nCRT cohorts)", n_primary),
                pool_override = pool_from_disk("IFN"))
# CTL 作支持性证据置左, NK(边界显著)置右, 不抢主线
pB1 <- forest_gg(mk_df("Cyto"),
                 title = sprintf("B  Cytotoxic lymphocytes (%d nCRT cohorts)", n_primary),
                 pool_override = pool_from_disk("Cyto"))
pB2 <- forest_gg(mk_df("NK"),
                 title = sprintf("   NK cells (%d nCRT cohorts)", n_primary),
                 show_ylab = FALSE, pool_override = pool_from_disk("NK"))
pB <- arrangeGrob(pB1, pB2, ncol = 2)

# ---------- C: IG 分组 CR 率 ----------
ig <- read.csv("GSE209746_ig_groups.csv", stringsAsFactors = FALSE)
ig$cr <- ig$response == "good"
tab <- table(ig$ig_group)
cr <- tapply(ig$cr, ig$ig_group, mean)
ncr <- tapply(ig$cr, ig$ig_group, sum)
c_df <- data.frame(IG = names(cr), n = as.integer(tab),
                   cr_n = as.integer(ncr), rate = as.numeric(cr))
c_df$IG <- factor(c_df$IG, levels = c("IG1", "IG2", "IG3", "IG4"))
hot <- ig$ig_group %in% c("IG3", "IG4")
p_ig <- fisher.test(table(ifelse(hot, "IG3/4", "IG1/2"), ig$cr))$p.value

pC <- ggplot(c_df, aes(x = IG, y = rate, fill = IG)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d/%d\n(%.0f%%)", cr_n, n, rate * 100),
                y = rate + 0.07), size = 3.2, family = TXT) +
  annotate("segment", x = 1.5, xend = 1.5, y = 0.80, yend = 0.88, linewidth = 0.4) +
  annotate("segment", x = 4.5, xend = 4.5, y = 0.80, yend = 0.88, linewidth = 0.4) +
  annotate("segment", x = 1.5, xend = 4.5, y = 0.88, yend = 0.88, linewidth = 0.4) +
  annotate("text", x = 3, y = 0.96,
           label = sprintf("IG3/4 vs IG1/2: Fisher P = %.3f", p_ig),
           size = 3.2, family = TXT) +
  scale_fill_manual(values = c(IG1 = "#9ECAE1", IG2 = "#4292C6",
                               IG3 = "#EF6548", IG4 = "#7F0000")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05), expand = c(0, 0)) +
  labs(title = "C  Sustained CR rate by immune group\n     (GSE209746)",
       x = "Immune group (IG1 cold to IG4 hot/MSI)", y = "2-year sustained CR rate") +
  theme_classic(base_size = 11, base_family = TXT) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        legend.position = "none",
        axis.text = element_text(color = "black"))

# ---------- D: GSE87211 DFS by pCR ----------
eset3 <- GEOquery:::parseGSEMatrix("GSE87211_series_matrix.txt.gz",
                                   AnnotGPL = FALSE, getGPL = FALSE)$eset
clin3 <- pData(eset3)
tissue_col <- grep("^tissue", colnames(clin3), value = TRUE)[1]
ypT_col <- grep("depth of invasion after", colnames(clin3), value = TRUE)[1]
ypN_col <- grep("lymph node metastasis after", colnames(clin3), value = TRUE)[1]
is_tum <- clin3[[tissue_col]] == "rectal tumor"
num_clin <- function(x) suppressWarnings(as.numeric(as.character(x)))
dfs_t <- num_clin(clin3[["disease free time (month):ch1"]])
dfs_e <- num_clin(clin3[["cancer recurrance after surgery:ch1"]])
ypT <- suppressWarnings(as.numeric(as.character(clin3[[ypT_col]])))
ypN <- suppressWarnings(as.numeric(as.character(clin3[[ypN_col]])))
pcr_all <- ifelse(ypT == 0 & ypN == 0, "pCR (ypT0N0)",
                  ifelse((!is.na(ypT) & ypT != 0) | (!is.na(ypN) & ypN != 0),
                         "non-pCR", NA))
ok <- is_tum & !is.na(dfs_t) & !is.na(dfs_e) & dfs_t > 0 & !is.na(pcr_all)
srvd <- data.frame(time = dfs_t[ok], event = dfs_e[ok], pCR = pcr_all[ok])
km <- survfit(Surv(time, event) ~ pCR, data = srvd)
sdf <- survdiff(Surv(time, event) ~ pCR, data = srvd)
p_km <- pchisq(sdf$chisq, df = 1, lower.tail = FALSE)

strata_vec <- rep(names(km$strata), km$strata)   # 每行所属分层
km_df <- data.frame(arm = sub("^pCR=", "", strata_vec),
                    time = km$time, surv = km$surv, nc = km$n.censor)
cens_df <- km_df[km_df$nc > 0, , drop = FALSE]
km_df <- rbind(km_df[, c("arm", "time", "surv")],
               data.frame(arm = unique(km_df$arm), time = 0, surv = 1))
lv <- c("pCR (ypT0N0)", "non-pCR")
km_df$arm <- factor(km_df$arm, levels = lv)
cens_df$arm <- factor(cens_df$arm, levels = lv)
ev <- table(srvd$pCR, factor(srvd$event, levels = c(0, 1)))

pD <- ggplot(km_df, aes(x = time, y = surv, color = arm)) +
  geom_step(linewidth = 0.9) +
  geom_point(data = cens_df, aes(x = time, y = surv),
             shape = 3, size = 1.5, stroke = 0.6, show.legend = FALSE) +
  annotate("text", x = 5, y = 0.33, hjust = 0, size = 3.2, family = TXT,
           label = sprintf("log-rank P = %.3f\npCR recurrences: %d/%d\nnon-pCR: %d/%d",
                           p_km, ev["pCR (ypT0N0)", "1"], sum(ev["pCR (ypT0N0)", ]),
                           ev["non-pCR", "1"], sum(ev["non-pCR", ]))) +
  scale_color_manual(values = c("pCR (ypT0N0)" = "#0072B2", "non-pCR" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
  labs(title = "D  DFS by pathological response (GSE87211)",
       x = "Months after surgery", y = "Disease-free survival", color = NULL) +
  theme_classic(base_size = 11, base_family = TXT) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        legend.position = c(0.60, 0.16),
        legend.text = element_text(size = 9.5),
        axis.text = element_text(color = "black"))

# ---------- 输出 ----------
ggsave("fig_A_forest_IFN_all.png", pA, width = 9.4, height = 5.6, dpi = 300)
ggsave("fig_B_forest_immune_all.png", pB, width = 13.6, height = 5.8, dpi = 300)
ggsave("fig_C_ig_crrate.png", pC, width = 5.4, height = 4.6, dpi = 300)
ggsave("fig_D_km_pcr.png", pD, width = 5.6, height = 4.6, dpi = 300)

# ---------- 组合主图: A 顶行 / B 中行 / C+D 底行 ----------
main_grob <- arrangeGrob(pA, pB, arrangeGrob(pC, pD, ncol = 2),
                        ncol = 1, heights = c(1.25, 1.15, 1.0))
png("Fig1_main.png", width = 3800, height = 4400, res = 300)
grid.draw(main_grob)
dev.off()

# 可编辑矢量输出，约为期刊双栏宽度。
svg("Fig1_main.svg", width = 7.2, height = 8.3, family = TXT)
grid.draw(main_grob)
dev.off()
cairo_pdf("Fig1_main.pdf", width = 7.2, height = 8.3, family = TXT)
grid.draw(main_grob)
dev.off()

cat("===== 完成 =====\n")
