# ==============================
# run_all.R — 全流程复现主脚本
# 项目: LARC nCRT 疗效的九队列主分析 + 放疗/化疗方案敏感性分析
# 用法: Rscript run_all.R [--preflight-only]
#       (可从任意工作目录调用, 原始数据须已就位)
#
# 阶段 0 (数据获取) 不在本脚本内:
#   - GEO series matrix / 补充计数表经 curl/GEOquery 下载 (README.md 第 3 节)
#   - 新队列原始文件置于 new_cohorts/
# 调试脚本 (debug_*.R) 不属于正式流程, 不纳入。
# ==============================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
project_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(project_dir)
source("project_config.R")
task_args <- commandArgs(trailingOnly = TRUE)
preflight_only <- "--preflight-only" %in% task_args

run_log <- file.path(project_dir, paste0("run_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_connection <- file(run_log, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message", append = TRUE)
warnings_seen <- character()
cat("Project directory: ", project_dir, "\n", sep = "")
cat("Run started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Configuration version: ", LARC_CONFIG$version, "\n", sep = "")

run_step <- function(f) {
  cat("\n", strrep("=", 60), "\n>>> ", f, "\n", strrep("=", 60), "\n", sep = "")
  t0 <- Sys.time()
  # Evaluate each step in its own child environment. Steps communicate through
  # files, so isolating objects prevents a previous script's variables from
  # masking missing-input or naming bugs in a later script.
  step_env <- new.env(parent = globalenv())
  # Explicit encoding conversion can silently return zero expressions when R
  # starts in the C locale on Windows. Parse first with the native UTF-8-aware
  # reader, reject empty parses, then evaluate in the isolated child environment.
  expr <- parse(file = f, keep.source = FALSE)
  if (!length(expr)) stop("No executable expressions parsed from ", f)
  withCallingHandlers(
    for (i in seq_along(expr)) eval(expr[[i]], envir = step_env),
    warning = function(w) {
      warning_text <- paste0(f, ": ", conditionMessage(w))
      warnings_seen <<- c(warnings_seen, warning_text)
      message("WARNING [", warning_text, "]")
      invokeRestart("muffleWarning")
    }
  )
  cat(sprintf("<<< %s 完成 (%.1f min)\n", f, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

pipeline <- c(
  # ---- 阶段 1: 发现队列 GSE209746 ----
  "GSE209746_alignment_check.R",     # 94→93 样本核对
  "GSE209746_DESeq2_analysis.R",     # DESeq2 + 15基因签名
  "GSE209746_loose_filter.R",        # 宽松过滤版 (保留低表达签名基因)
  "GSE209746_MSS_only.R",            # MSS-only 敏感性
  "GSE209746_immune_groups.R",       # IG1-4 免疫分组重建
  "IFN_signature_ROC_MSI.R",         # IFN 分数 ROC + MSI 分析
  "GSE209746_genomic_adjust.R",      # WES 纯度/TMB/CNA/GD 校正

  # ---- 阶段 2: 既有验证队列 ----
  "GSE35452_analysis.R",
  "GSE35452_crosscheck.R",
  "GSE87211_analysis.R",
  "GSE150082_validation.R",
  "GSE45404_analysis.R",
  "hallmark_3cohorts.R",             # camera Hallmark 一致性
  "immune_ssgsea_3cohorts.R",        # ssGSEA 免疫打分 (与 Bindea 对标)
  "tumor_purity_adjust.R",           # MCP 间质校正
  "combined_model.R",                # IFN+IGF2+L1CAM 固定系数转移模型
  "GSE87211_survival.R",
  "GSE87211_survival_v2.R",          # DFS/OS + yp 校正 + 样条 (正式版)
  "TCGA_READ_survival.R",

  # ---- 阶段 3: 系统检索新增队列 (原始数据在 new_cohorts/) ----
  "check_eligibility.R",             # 41 候选资格评估 → eligibility_report.txt
  "new_cohorts_process1.R",          # GSE119409 / GSE213331
  "new_cohorts_process3.R",          # GSE53781 + GSE94104 + GSE133057 注释回填
  "new_cohorts_process4.R",          # GSE56699：重复活检按患者合并（放疗敏感性）
  "validate_outputs.R",              # 队列数、疗效标签和核心特征防呆校验

  # ---- 阶段 4: 九队列 nCRT 主 meta + 非 nCRT 方案敏感性分析 ----
  "meta_all_cohorts.R",              # nCRT 主分析; 放疗/化疗/全部方案敏感性
  "meta_endpoint_subgroups.R",       # CR/pCR vs broader response-definition 亚组
  "meta_diagnostics.R",              # LOO + Egger + 漏斗图
  "meta_baujat.R",                   # Baujat 异质性贡献与影响诊断
  "immune_hot_xcohort.R",            # 队列内 top-tertile 分数敏感性 OR 合并
  "GSE87211_clinical_utility.R",     # cT/cN+IFN 临床效用 + KRAS + 反向KM随访
  "gse87211_regimen_sensitivity.R",  # 队列内 5-FU vs ±Oxa 方案敏感性
  "genes15_meta.R",                  # 15 基因逐基因 meta
  "meta_prediction_intervals.R",     # 95% 预测区间 (IntHout 2016)
  "ifn_robustness.R",                # 替代 IFN 签名与核心单基因稳健性

  # ---- 阶段 5: 图 ----
  "figures_main.R",                  # Fig1 + 面板 A-D
  "fig_S0_prisma.R",                 # PRISMA 流程图
  "fig_S2_diagnostics.R",            # LOO + 漏斗图 + Baujat 三行组版
  "fig_supplement_assemblies.R",      # Figure S4 / S6 双面板组版
  "assemble_supplementary_tables.R"   # 投稿版 Table S6 / S7
)

# Fail early with actionable diagnostics instead of producing a partial set of
# stale outputs after a missing package or source file is encountered.
required_packages <- c(
  "GEOquery", "GSVA", "Biobase", "DESeq2", "limma", "pROC", "survival",
  "ggplot2", "gridExtra", "org.Hs.eg.db", "illuminaHumanv3.db", "scales",
  "svglite", "png"
)
required_inputs <- c(
  "h.all.symbols.gmt", "reactome2022.gmt", "mcpcounter_genes.txt",
  "GSE209746_geneExp_n_114_ensembl.csv.gz", "GSE209746_series_matrix.txt",
  "paper_94_labels.csv", "paper_S8A_results.csv", "paper_S1C_wes.csv",
  "paper_S1D_bindea_scores.csv", "GSE35452_eset.rds",
  "GSE87211_series_matrix.txt.gz", "GSE150082_series_matrix.txt.gz",
  "GSE45404-GPL570_series_matrix.txt.gz", "GPL13497_probe_map.tsv",
  "TCGA-READ.star_counts.tsv.gz", "READ_survival.txt", "READ_clinicalMatrix",
  "new_cohorts/GPL6102.annot.gz",
  "new_cohorts/GPL18134.soft.gz", "new_cohorts/GSE119409_series_matrix.txt.gz",
  "new_cohorts/GSE133057_series_matrix.txt.gz", "new_cohorts/GSE213331_series_matrix.txt.gz",
  "new_cohorts/GSE213331_readcount_matrix.txt.gz", "new_cohorts/GSE53781_series_matrix.txt.gz",
  "new_cohorts/GSE94104_series_matrix.txt.gz",
  "new_cohorts/E-GEOD-56699.processed.1.zip",
  "new_cohorts/E-GEOD-56699.sdrf.txt",
  "new_cohorts/A-GEOD-14951.adf.txt"
)
required_support_scripts <- c("project_config.R", "meta_utils.R", "cohort_utils.R")

preflight <- function() {
  missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  missing_scripts <- pipeline[!file.exists(pipeline)]
  missing_support <- required_support_scripts[!file.exists(required_support_scripts)]
  if (length(missing_packages)) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
         ". Install them before running the pipeline.")
  }
  if (length(missing_inputs)) {
    stop("Missing required input files: ", paste(missing_inputs, collapse = ", "))
  }
  if (length(missing_scripts)) {
    stop("Missing pipeline scripts: ", paste(missing_scripts, collapse = ", "))
  }
  if (length(missing_support)) {
    stop("Missing support scripts: ", paste(missing_support, collapse = ", "))
  }
  cat("Preflight PASS: ", length(required_packages), " packages, ",
      length(required_inputs), " input files, ", length(pipeline),
      " pipeline scripts and ", length(required_support_scripts),
      " support scripts are available.\n", sep = "")
  invisible(TRUE)
}

pipeline_ok <- FALSE
pipeline_error <- NULL
tryCatch({
  preflight()
  if (preflight_only) {
    cat("--preflight-only specified; analysis steps were not executed.\n")
  } else {
    for (f in pipeline) run_step(f)
    session_header <- trimws(c(
      "Session info for: nine-cohort nCRT primary meta-analysis with two non-nCRT sensitivity cohorts",
      paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "",
      capture.output(sessionInfo())
    ), which = "right")
    writeLines(session_header, "session_info.txt", useBytes = TRUE)
  }
  pipeline_ok <- TRUE
  if (preflight_only) {
    cat("\n===== 预检完成（未执行分析） =====\n")
  } else {
    cat("Pipeline warning count: ", length(warnings_seen), "\n", sep = "")
    cat("\n===== 全流程完成 =====\n")
  }
}, error = function(e) {
  pipeline_error <<- conditionMessage(e)
  message("Pipeline failed: ", pipeline_error)
}, finally = {
  cat("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  sink(type = "message")
  sink()
  close(log_connection)
})

if (!pipeline_ok) stop("Pipeline failed; see ", run_log, ": ", pipeline_error, call. = FALSE)
cat("Run log: ", run_log, "\n", sep = "")
