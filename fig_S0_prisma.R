# ==============================
# Figure S0: PRISMA-style 队列检索筛选流程图
# 数字来源: cohort_search_record.md (2026-09-01 检索；2026-09-10 检索后复核)
# 输出: fig_S0_prisma.png / fig_S0_prisma.svg
# ==============================

suppressMessages(library(grid))

box_grob <- function(x, y, w, h, label, fill = "grey97", fs = 9.5, bold_first = FALSE) {
  grobTree(
    roundrectGrob(x = x, y = y, width = w, height = h,
                  r = unit(0.012, "snpc"),
                  gp = gpar(fill = fill, col = "grey35", lwd = 0.8)),
    textGrob(label, x = x, y = y, gp = gpar(fontsize = fs), just = "center")
  )
}
arrow_down <- function(x, y0, y1) {
  segmentsGrob(x, y0, x, y1, gp = gpar(lwd = 1.1),
               arrow = arrow(type = "closed", length = unit(0.09, "inches")))
}
arrow_right <- function(x0, x1, y) {
  segmentsGrob(x0, y, x1, y, gp = gpar(lwd = 1.1),
               arrow = arrow(type = "closed", length = unit(0.09, "inches")))
}

png("fig_S0_prisma.png", width = 3000, height = 2100, res = 300)
grid.newpage()

# 左侧竖排三框 (identification -> assessed -> included)
bx <- 0.30; bw <- 0.44
b1y <- 0.86; b2y <- 0.56; b3y <- 0.18; bh <- 0.155

grid.draw(box_grob(bx, b1y, bw, bh,
  "Records identified\nGEO DataSets search (n = 38)\nInitial citation chasing (n = 2)\nPost-search audit (n = 1)",
  fill = "#DEEBF7"))
grid.draw(box_grob(bx, b2y, bw, bh,
  "Datasets assessed for eligibility\n(n = 41)"))
grid.draw(box_grob(bx, b3y, bw, bh,
  "Cohorts included in quantitative synthesis\n9 nCRT cohorts (n = 577)\n+ 1 radiotherapy sensitivity cohort (n = 51)",
  fill = "#FDE9D9"))

grid.draw(arrow_down(bx, b1y - bh/2, b2y + bh/2))
grid.draw(arrow_down(bx, b2y - bh/2, b3y + bh/2))

# 右侧排除框
ex <- 0.78; ew <- 0.38
grid.draw(box_grob(ex, (b2y + b3y)/2, ew, 0.24,
  "Excluded (n = 31)\nResponse labels not deposited (n = 4)\nPlatform annotation unavailable (n = 2)\nDuplicate cohort (n = 1)\nIneligible design (n = 24)",
  fill = "grey95"))
grid.draw(arrow_right(bx + bw/2, ex - ew/2, b2y - bh/2 - 0.035))

# 阶段标签
grid.draw(textGrob("Identification", x = 0.045, y = b1y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))
grid.draw(textGrob("Screening", x = 0.045, y = b2y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))
grid.draw(textGrob("Included", x = 0.045, y = b3y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))

# 脚注
grid.draw(textGrob("Search date: 2026-09-01; post-search audit: 2026-09-10.  GSE56699 radiotherapy-only. GSE213331 excluded (post-treatment resections).\nTCGA-READ (n = 164) is prognostic context only.",
                   x = 0.5, y = 0.045, gp = gpar(fontsize = 8, col = "grey35")))

dev.off()

svg("fig_S0_prisma.svg", width = 10, height = 7)
grid.newpage()
grid.draw(box_grob(bx, b1y, bw, bh,
  "Records identified\nGEO DataSets search (n = 38)\nInitial citation chasing (n = 2)\nPost-search audit (n = 1)",
  fill = "#DEEBF7"))
grid.draw(box_grob(bx, b2y, bw, bh,
  "Datasets assessed for eligibility\n(n = 41)"))
grid.draw(box_grob(bx, b3y, bw, bh,
  "Cohorts included in quantitative synthesis\n9 nCRT cohorts (n = 577)\n+ 1 radiotherapy sensitivity cohort (n = 51)",
  fill = "#FDE9D9"))
grid.draw(arrow_down(bx, b1y - bh/2, b2y + bh/2))
grid.draw(arrow_down(bx, b2y - bh/2, b3y + bh/2))
grid.draw(box_grob(ex, (b2y + b3y)/2, ew, 0.24,
  "Excluded (n = 31)\nResponse labels not deposited (n = 4)\nPlatform annotation unavailable (n = 2)\nDuplicate cohort (n = 1)\nIneligible design (n = 24)",
  fill = "grey95"))
grid.draw(arrow_right(bx + bw/2, ex - ew/2, b2y - bh/2 - 0.035))
grid.draw(textGrob("Identification", x = 0.045, y = b1y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))
grid.draw(textGrob("Screening", x = 0.045, y = b2y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))
grid.draw(textGrob("Included", x = 0.045, y = b3y, rot = 90,
                   gp = gpar(fontsize = 10, fontface = "bold", col = "grey40")))
grid.draw(textGrob("Search date: 2026-09-01; post-search audit: 2026-09-10.  GSE56699 radiotherapy-only. GSE213331 excluded (post-treatment resections).\nTCGA-READ (n = 164) is prognostic context only.",
                   x = 0.5, y = 0.045, gp = gpar(fontsize = 8, col = "grey35")))
dev.off()
cat("输出: fig_S0_prisma.png / fig_S0_prisma.svg\n")
