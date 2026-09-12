# =====================================================================
# Assemble raster source plots into consistently labelled supplementary panels.
# Source plots are regenerated earlier in run_all.R; this step only composes them.
# Outputs: fig_S4_models.* and fig_S6_survival.*
# =====================================================================

suppressMessages({
  library(grid)
  library(png)
})

panel_grob <- function(path, label) {
  if (!file.exists(path)) stop("Missing panel source: ", path)
  grobTree(
    rasterGrob(readPNG(path), width = unit(1, "npc"), height = unit(1, "npc"),
               interpolate = TRUE),
    textGrob(label, x = unit(0.012, "npc"), y = unit(0.985, "npc"),
             just = c("left", "top"),
             gp = gpar(fontface = "bold", fontsize = 13,
                       col = "black", fill = "white"))
  )
}

make_two_panel <- function(files, stem) {
  grob <- gTree(children = gList(
    gTree(children = gList(panel_grob(files[[1]], "A")),
          vp = viewport(x = 0.25, width = 0.5)),
    gTree(children = gList(panel_grob(files[[2]], "B")),
          vp = viewport(x = 0.75, width = 0.5))
  ))
  width_in <- 183 / 25.4
  height_in <- 92 / 25.4
  png(paste0(stem, ".png"), width = width_in, height = height_in,
      units = "in", res = 300, bg = "white")
  grid.newpage(); grid.draw(grob); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in,
                       height = height_in, family = "Arial")
  grid.newpage(); grid.draw(grob); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  grid.newpage(); grid.draw(grob); dev.off()
}

make_two_panel(
  c("IFN_score_ROC_3cohorts.png", "combined_model_ROC.png"),
  "fig_S4_models"
)
make_two_panel(
  c("GSE87211_IFN_DFS_spline.png", "TCGA_READ_IFN_OS_KM.png"),
  "fig_S6_survival"
)

cat("Outputs: fig_S4_models.png/.pdf/.svg; fig_S6_survival.png/.pdf/.svg\n")
cat("===== Complete =====\n")
