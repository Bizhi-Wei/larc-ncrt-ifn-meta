# ==============================
# GSE209746 免疫分组(IG1-IG4)复现 + MSS-only 敏感性分析
# 方法(按论文 Nature Med PMID 35970919): IG4 = 全部 MSI 单独成组;
#   MSS 样本用论文 Table S1D 的 22 个 Bindea 免疫 ssGSEA 分数 kmeans(k=3),
#   按免疫热度排序 IG1(冷)/IG2(中)/IG3(热)
# ==============================

suppressMessages({ library(DESeq2) })

vsd <- readRDS("GSE209746_vsd.rds")
labels <- read.csv("paper_94_labels.csv", stringsAsFactors = FALSE)
s1d <- read.csv("paper_S1D_bindea_scores.csv", stringsAsFactors = FALSE)
sams <- colnames(assay(vsd))
grp <- ifelse(labels$response[match(sams, labels$patient)] == "CR", "good", "poor")
names(grp) <- sams
msi <- ifelse(setNames(s1d$MSI_status, s1d$patient)[sams] == 1, "MSI", "MSS")
msi[is.na(msi)] <- "NA"; names(msi) <- sams

bindea_imm <- c("B_CELLS","CD8_T_CELLS","CYTOTOXIC_CELLS","DC","EOSINOPHILS","IDC",
                "MACROPHAGES","MAST_CELLS","NEUTROPHILS","NK_CD56BRIGHT_CELLS",
                "NK_CD56DIM_CELLS","NK_CELLS","TCM","TEM","TFH","TGD","TH17_CELLS",
                "TH1_CELLS","TH2_CELLS","TREG","T_CELLS","T_HELPER_CELLS")
s1d_m <- s1d[match(sams, s1d$patient), ]
mm <- sapply(s1d_m[, bindea_imm], as.numeric)   # 93样本 x 22细胞类型
z <- t(scale(mm))                                # 按细胞类型z化后转置: 22 x 93
colnames(z) <- sams

mss_sams <- sams[msi == "MSS"]
zm <- z[, mss_sams]

# ---- 种子稳定性扫描: kmeans(k=3) 在 50 个种子下的 hot/cold CR 率与 P ----
set_list <- list()
for (sd in 1:50) {
  set.seed(sd)
  cl <- kmeans(t(zm), centers = 3, nstart = 50)$cluster
  ord <- order(sapply(split(colMeans(zm), cl), mean))
  ig <- paste0("IG", ord[cl]); names(ig) <- mss_sams
  g <- grp[mss_sams]
  hot <- ig == "IG3"; cold <- ig %in% c("IG1", "IG2")
  set_list[[sd]] <- data.frame(seed = sd, n_ig3 = sum(hot),
                               cr_ig3 = mean(g[hot] == "good"),
                               cr_cold = mean(g[cold] == "good"),
                               p = fisher.test(table(ifelse(hot, "hot", "cold"), g))$p.value)
}
stab <- do.call(rbind, set_list)
cat("===== kmeans 种子稳定性(50 个种子, MSS-only IG3 vs IG1/2) =====\n")
cat(sprintf("IG3 大小: 中位 %d (范围 %d-%d)\n", median(stab$n_ig3), min(stab$n_ig3), max(stab$n_ig3)))
cat(sprintf("IG3 CR率: 中位 %.0f%% (范围 %.0f%%-%.0f%%)\n",
            median(stab$cr_ig3) * 100, min(stab$cr_ig3) * 100, max(stab$cr_ig3) * 100))
cat(sprintf("IG1/2 CR率: 中位 %.0f%% (范围 %.0f%%-%.0f%%)\n",
            median(stab$cr_cold) * 100, min(stab$cr_cold) * 100, max(stab$cr_cold) * 100))
cat(sprintf("Fisher P: 中位 %.4f; P<0.05 的种子占 %d/50\n", median(stab$p), sum(stab$p < 0.05)))

# ---- 代表解: seed=42(与复现一致的那个) ----
set.seed(42)
cl <- kmeans(t(zm), centers = 3, nstart = 50)$cluster
ord <- order(sapply(split(colMeans(zm), cl), mean))
ig_mss <- paste0("IG", ord[cl]); names(ig_mss) <- mss_sams
ig_all <- c(ig_mss,
            setNames(rep("IG4", sum(msi == "MSI")), sams[msi == "MSI"]))
ig_all <- ig_all[sams[!is.na(grp)]]
ig_all <- ig_all[!is.na(ig_all)]   # 3 例 MSI 状态缺失样本不进入分组分析
cat("进入分组分析样本:", length(ig_all), "(85 MSS + 5 MSI; 3 例 MSI-NA 剔除)\n")
g <- grp[names(ig_all)]

cat("\n===== 全队列(IG3/IG4 vs IG1/IG2, 含MSI) =====\n")
hot <- ig_all %in% c("IG3", "IG4"); cold <- ig_all %in% c("IG1", "IG2")
tab_full <- table(group = ifelse(hot, "IG3/4", "IG1/2"), response = g)
print(table(IG = ig_all, response = g))
print(tab_full)
ft <- fisher.test(tab_full)
cat(sprintf("CR率 IG3/4 = %.0f%% vs IG1/2 = %.0f%%, Fisher P = %.4f\n",
            mean(g[hot] == "good") * 100, mean(g[cold] == "good") * 100, ft$p.value))

cat("\n===== MSS-only(IG3 vs IG1/IG2, 剔除IG4/MSI) =====\n")
keep <- names(ig_all)[ig_all != "IG4"]
ig2 <- ig_all[keep]; g2 <- grp[keep]
hot2 <- ig2 == "IG3"; cold2 <- ig2 %in% c("IG1", "IG2")
tab_mss <- table(group = ifelse(hot2, "IG3", "IG1/2"), response = g2)
print(tab_mss)
ft2 <- fisher.test(tab_mss)
cat(sprintf("CR率 IG3 = %.0f%% vs IG1/2 = %.0f%%, Fisher P = %.4f\n",
            mean(g2[hot2] == "good") * 100, mean(g2[cold2] == "good") * 100, ft2$p.value))

cat("\n===== IG4(MSI) 单独 =====\n")
g4 <- grp[sams[msi == "MSI"]]
g4 <- g4[!is.na(g4)]
cat(sprintf("IG4(MSI) n=%d, CR=%d (%.0f%%)\n", length(g4), sum(g4 == "good"),
            mean(g4 == "good") * 100))

# ---- 保存分组与汇总 ----
write.csv(data.frame(patient = names(ig_all), ig_group = ig_all,
                     response = grp[names(ig_all)], msi = msi[names(ig_all)]),
          "GSE209746_ig_groups.csv", row.names = FALSE)
summ <- data.frame(
  analysis = c("full_IG34_vs_IG12", "MSS_only_IG3_vs_IG12", "IG4_MSI_only"),
  n_hot = c(sum(hot), sum(hot2), length(g4)),
  cr_hot = c(mean(g[hot] == "good"), mean(g2[hot2] == "good"), mean(g4 == "good")),
  n_cold = c(sum(cold), sum(cold2), NA),
  cr_cold = c(mean(g[cold] == "good"), mean(g2[cold2] == "good"), NA),
  fisher_p = c(ft$p.value, ft2$p.value, NA))
write.csv(summ, "GSE209746_ig_groups_summary.csv", row.names = FALSE)

cat("\n===== 完成 =====\n")
