# ==============================
# 候选队列资格筛查: 解析各 series matrix 的表型字段
# ==============================
suppressMessages(library(GEOquery))

gses <- c("GSE216616","GSE139255","GSE40492","GSE46862","GSE119409","GSE318268",
          "GSE133057","GSE213331","GSE53781","GSE242786","GSE93375","GSE60331","GSE15781")

report_lines <- capture.output({
for (g in gses) {
  f <- paste0("new_cohorts/", g, "_series_matrix.txt.gz")
  cat("\n", rep("=", 20), g, rep("=", 20), "\n", sep = "")
  eset <- tryCatch(
    GEOquery:::parseGSEMatrix(f, AnnotGPL = FALSE, getGPL = FALSE)$eset,
    error = function(e) { cat("PARSE ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(eset)) next
  pd <- pData(eset)
  cat("samples:", nrow(pd), " features:", nrow(exprs(eset)), "\n")
  ch <- grep("characteristics_ch1", colnames(pd), value = TRUE)
  cat("-- characteristics 字段 --\n")
  for (cc in ch) {
    vals <- as.character(pd[[cc]])
    key <- sub(":.*", "", vals[!is.na(vals)][1])
    uv <- unique(sub("^[^:]*: *", "", vals))
    if (length(uv) > 12) uv <- c(head(uv, 6), sprintf("... (%d unique)", length(uv)))
    cat(sprintf("  %s [%s]: %s\n", cc, key, paste(uv, collapse = " | ")))
  }
  for (extra in c("source_name_ch1", "title", "treatment_protocol_ch1")) {
    if (extra %in% colnames(pd)) {
      uv <- unique(as.character(pd[[extra]]))
      if (length(uv) > 8) uv <- c(head(uv, 4), sprintf("... (%d unique)", length(uv)))
      cat(sprintf("  %s: %s\n", extra, paste(uv, collapse = " | ")))
    }
  }
}
})
writeLines(sub("[ \\t]+$", "", report_lines),
           "new_cohorts/eligibility_report.txt", useBytes = TRUE)
cat("完成: new_cohorts/eligibility_report.txt\n")
