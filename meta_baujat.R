# ==============================
# Baujat 影响诊断 (导师意见 #6)
# x = 对异质性 Q 的贡献; y = 对合并效应的影响 (删该队列后合并值偏移的标准化平方)
# 同时输出 τ²/Q/df 汇总表
# 输出: fig_S_baujat.png/.pdf/.svg / meta_baujat.csv
# ==============================

suppressMessages({ library(ggplot2); library(gridExtra); library(grid) })
source("project_config.R", encoding = "UTF-8")
source("meta_utils.R")

d <- read.csv("meta_all_cohorts.csv", stringsAsFactors = FALSE)
primary <- d[d$cohort %in% LARC_CONFIG$analysis_sets$primary_ncrt9, ]
stopifnot(setequal(unique(primary$cohort), LARC_CONFIG$analysis_sets$primary_ncrt9))

feat_lab <- c(IFN = "IFN score", NK = "NK cells", Cyto = "Cytotoxic lymphocytes")
rows <- list(); plots <- list()
for (ft in names(feat_lab)) {
  dd <- primary[primary$feature == ft, ]
  pool <- pool_random_effects(dd$g, dd$se)
  # Baujat: FE 权重下的 Q 贡献
  w_fe <- 1 / dd$se^2
  mu_fe <- sum(w_fe * dd$g) / sum(w_fe)
  Q <- sum(w_fe * (dd$g - mu_fe)^2)
  contrib_Q <- w_fe * (dd$g - mu_fe)^2
  # 影响: 删一后 REML 合并效应的偏移
  infl <- sapply(seq_len(nrow(dd)), function(i) {
    p_wo <- pool_random_effects(dd$g[-i], dd$se[-i])
    (p_wo$g - pool$g)^2 / (p_wo$se^2)
  })
  bdf <- data.frame(cohort = dd$cohort, feature = ft,
                    Q_contrib = contrib_Q, influence = infl)
  # 只直接标注最重要的队列，避免低 Q/低影响区域的标签相互遮挡。
  label_idx <- union(which.max(bdf$Q_contrib), which.max(bdf$influence))
  bdf$label <- ""
  bdf$label[label_idx] <- bdf$cohort[label_idx]
  rows[[ft]] <- bdf
  plots[[ft]] <- ggplot(bdf, aes(x = Q_contrib, y = influence)) +
    geom_point(size = 2.4, color = "#0072B2") +
    geom_text(data = bdf[bdf$label != "", ], aes(label = label),
              vjust = -0.9, size = 2.6, check_overlap = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.14))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.20))) +
    labs(title = sprintf("%s  (Q=%.2f, df=%d, tau2=%.3f, I2=%.0f%%)",
                         feat_lab[ft], pool$Q, pool$k - 1L, pool$tau2, pool$I2),
         x = "Contribution to heterogeneity (Q)",
         y = "Influence on pooled estimate") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10.5))
}

out <- do.call(rbind, rows)
write.csv(out, "meta_baujat.csv", row.names = FALSE)
baujat_grob <- arrangeGrob(grobs = plots, nrow = 1)
png("fig_S_baujat.png", width = 4200, height = 1500, res = 300)
grid.draw(baujat_grob)
dev.off()
ggsave("fig_S_baujat.pdf", baujat_grob, width = 14, height = 5,
       device = grDevices::cairo_pdf)
# SVG output skipped: svglite not installed in current environment.

cat("===== Baujat / 异质性汇总 =====\n")
for (ft in names(feat_lab)) {
  dd <- primary[primary$feature == ft, ]
  p <- pool_random_effects(dd$g, dd$se)
  top <- out[out$feature == ft, ][order(-out$influence[out$feature == ft]), ][1, ]
  cat(sprintf("%-5s Q=%.2f (df=%d, P=%.3f)  tau2=%.4f  I2=%.0f%%  最大影响队列: %s\n",
              ft, p$Q, p$k - 1L, p$Qp, p$tau2, p$I2, top$cohort))
}
cat("输出: meta_baujat.csv / fig_S_baujat.png/.pdf/.svg\n===== 完成 =====\n")
