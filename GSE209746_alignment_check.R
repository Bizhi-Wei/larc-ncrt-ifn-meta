# ==============================
# GSE209746 表达矩阵与临床标签对齐验证
# (不做差异分析)
# ==============================

library(GEOquery)

# ---- 1. 读取表达矩阵 ----
exp_raw <- as.matrix(
  read.csv(
    gzfile("GSE209746_geneExp_n_114_ensembl.csv.gz"),
    row.names = 1,
    check.names = FALSE
  )
)

cat("===== dim(exp_raw) =====\n")
print(dim(exp_raw))

cat("\n===== colnames(exp_raw)[1:10] =====\n")
print(colnames(exp_raw)[1:10])

# ---- 2. 获取临床信息(优先在线,失败则用本地文件) ----
gse <- tryCatch(
  getGEO("GSE209746", GSEMatrix = TRUE, AnnotGPL = FALSE),
  error = function(e) {
    cat("在线 getGEO 失败,改用本地 series matrix 文件\n")
    cat("错误信息:", conditionMessage(e), "\n")
    GEOquery:::parseGSEMatrix(
      "GSE209746_series_matrix.txt",
      AnnotGPL = FALSE,
      getGPL = FALSE
    )$eset
  }
)
eset <- if (is.list(gse)) gse[[1]] else gse
clinical <- pData(eset)

cat("\n===== dim(clinical) =====\n")
print(dim(clinical))

# ---- 3. 找 response 列并整理成 CR / iCR / NA ----
response_col <- grep("response", colnames(clinical),
                     ignore.case = TRUE, value = TRUE)
cat("\n===== 含 response 的列名 =====\n")
print(response_col)

# 逐列查看内容,挑出真正有 CR/iCR 标签的列
for (cl in response_col) {
  cat("\n--- 列:", cl, "---\n")
  print(table(clinical[[cl]], useNA = "ifany"))
}

# 用第一列做分组
response_raw <- as.character(clinical[[response_col[1]]])

response <- ifelse(
  grepl("\\biCR\\b", response_raw, ignore.case = TRUE),
  "iCR",
  ifelse(
    grepl("\\bCR\\b", response_raw, ignore.case = TRUE),
    "CR",
    NA
  )
)

cat("\n===== table(response, useNA='ifany') =====\n")
print(table(response, useNA = "ifany"))

# ---- 4. 检查表达矩阵列名与 clinical 的匹配 ----
exp_cols <- colnames(exp_raw)

n_match_title <- sum(exp_cols %in% clinical$title)
n_match_gsm   <- sum(exp_cols %in% clinical$geo_accession)

cat("\n===== 列名匹配情况 =====\n")
cat("表达矩阵总列数                :", length(exp_cols), "\n")
cat("能匹配 clinical$title 的列数   :", n_match_title, "\n")
cat("能匹配 clinical$geo_accession :", n_match_gsm, "\n")

# 若按 title 匹配,验证每个表达样本都能拿到标签且不错位
if (n_match_title > 0) {
  idx <- match(exp_cols, clinical$title)
  cat("未匹配到 title 的列:",
      paste(exp_cols[is.na(idx)], collapse = ", "), "\n")
  matched_response <- response[idx]
  cat("\n===== 按 title 匹配后,表达矩阵114列的标签分布 =====\n")
  print(table(matched_response, useNA = "ifany"))

  # 抽查前5列的对应关系
  cat("\n===== 抽查前5列的对应关系 =====\n")
  print(data.frame(
    exp_colname      = exp_cols[1:5],
    matched_title    = clinical$title[idx[1:5]],
    matched_gsm      = clinical$geo_accession[idx[1:5]],
    matched_response = matched_response[1:5]
  ))
}

cat("\n===== 完成 =====\n")
