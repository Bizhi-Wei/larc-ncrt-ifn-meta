# =====================================================================
# GSE56699: pretreatment rectal-cancer FFPE biopsies before radiotherapy
#
# ArrayExpress E-GEOD-56699 contains 58 pretreatment (PB) arrays, but five
# patients have repeated arrays (six extra arrays) and one unique patient has
# no Mandard grade.  The analysis unit is therefore the patient: repeated
# arrays are averaged on the log2 expression scale before ssGSEA, leaving 51
# independent patients with Mandard grades (1-2 = good; 3-5 = poor).
#
# This cohort is radiotherapy-only according to the source-study description;
# it is sensitivity-only and must not enter the nCRT primary analysis.
# =====================================================================

suppressMessages({
  library(GSVA)
  library(org.Hs.eg.db)
})
source("cohort_utils.R")

sig_sets <- load_signature_sets()
ifn_sets <- sig_sets$ifn_sets
mcp_sets <- sig_sets$mcp_sets

zip_path <- "new_cohorts/E-GEOD-56699.processed.1.zip"
sdrf_path <- "new_cohorts/E-GEOD-56699.sdrf.txt"
adf_path <- "new_cohorts/A-GEOD-14951.adf.txt"

sdrf <- read.delim(
  sdrf_path, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
required_sdrf <- c(
  "Assay Name", "Derived Array Data File",
  "Characteristics [sample type]", "Characteristics [mandard]",
  "Characteristics [response 3 classes]", "Characteristics [patient code]",
  "Characteristics [origin]", "Characteristics [diagnosis]"
)
stopifnot(nrow(sdrf) == 72L, all(required_sdrf %in% names(sdrf)))

sample_type <- trimws(sdrf[["Characteristics [sample type]"]])
mandard <- suppressWarnings(as.integer(sdrf[["Characteristics [mandard]"]]))
patient <- trimws(sdrf[["Characteristics [patient code]"]])
assay <- trimws(sdrf[["Assay Name"]])
data_file <- trimws(sdrf[["Derived Array Data File"]])
is_pre <- sample_type == "PB"
has_trg <- mandard %in% 1:5
keep <- is_pre & has_trg

stopifnot(
  sum(is_pre) == 58L,
  sum(keep) == 57L,
  length(unique(patient[is_pre])) == 52L,
  length(unique(patient[keep])) == 51L,
  !anyNA(patient[keep]), all(nzchar(patient[keep])),
  !anyDuplicated(assay[keep])
)

array_response <- ifelse(
  mandard %in% c(1L, 2L), "good",
  ifelse(mandard %in% 3:5, "poor", NA_character_)
)
mapping <- data.frame(
  assay = assay,
  patient_code = patient,
  sample_type = sample_type,
  mandard = mandard,
  response_3_classes = sdrf[["Characteristics [response 3 classes]"]],
  response_binary = array_response,
  included_array = keep,
  exclusion_reason = ifelse(
    !is_pre, "post-treatment surgical specimen",
    ifelse(!has_trg, "pretreatment biopsy without Mandard grade", "")
  ),
  origin = sdrf[["Characteristics [origin]"]],
  diagnosis = sdrf[["Characteristics [diagnosis]"]],
  stringsAsFactors = FALSE
)
mapping$patient_array_count <- ave(
  rep.int(1L, nrow(mapping)), mapping$patient_code, FUN = length
)
write.csv(mapping, "GSE56699_response_mapping.csv", row.names = FALSE)

zip_members <- unzip(zip_path, list = TRUE)$Name
stopifnot(all(data_file[keep] %in% zip_members))

read_processed <- function(member, expected_reporters = NULL) {
  con <- unz(zip_path, member, open = "rt")
  on.exit(close(con), add = TRUE)
  dat <- read.delim(
    con, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = ""
  )
  stopifnot(all(c("Reporter Identifier", "VALUE") %in% names(dat)))
  reporter <- as.character(dat[["Reporter Identifier"]])
  if (!is.null(expected_reporters)) stopifnot(identical(reporter, expected_reporters))
  value <- suppressWarnings(as.numeric(dat[["VALUE"]]))
  stopifnot(length(value) == length(reporter), all(is.finite(value)), all(value >= 0))
  list(reporter = reporter, value = value)
}

selected_files <- data_file[keep]
selected_assays <- assay[keep]
first <- read_processed(selected_files[1])
reporter <- first$reporter
y <- matrix(NA_real_, nrow = length(reporter), ncol = length(selected_files),
            dimnames = list(reporter, selected_assays))
y[, 1] <- first$value
if (length(selected_files) > 1L) {
  for (j in 2:length(selected_files)) {
    y[, j] <- read_processed(selected_files[j], reporter)$value
  }
}

# The deposited values are positive loess-normalized intensities.  The log2
# transform is monotone within each array (and therefore rank-preserving for
# ssGSEA) while making patient-level averaging of repeated arrays meaningful.
y <- log2(y + 1)

selected_patient <- patient[keep]
patient_index <- split(seq_along(selected_patient), selected_patient)
patient_y <- vapply(
  patient_index,
  function(idx) rowMeans(y[, idx, drop = FALSE]),
  numeric(nrow(y))
)
rownames(patient_y) <- rownames(y)

patient_trg <- vapply(
  patient_index,
  function(idx) {
    value <- unique(mandard[keep][idx])
    stopifnot(length(value) == 1L)
    value
  },
  integer(1)
)
grp <- ifelse(patient_trg %in% c(1L, 2L), "good", "poor")
names(grp) <- names(patient_trg)
stopifnot(
  ncol(patient_y) == 51L,
  sum(grp == "good") == 23L,
  sum(grp == "poor") == 28L
)

# A-GEOD-14951 supplies GenBank accessions for the DASL v4 probes.  Resolve
# those accessions through org.Hs.eg.db rather than relying on the invalid
# 990-byte GPL14951.annot.gz file that previously contained an HTML 404 page.
adf <- read.delim(
  adf_path, skip = 15L, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
required_adf <- c("Reporter Name", "Reporter Database Entry [genbank]")
stopifnot(all(required_adf %in% names(adf)), !anyDuplicated(adf[["Reporter Name"]]))
accession <- setNames(
  adf[["Reporter Database Entry [genbank]"]], adf[["Reporter Name"]]
)[rownames(patient_y)]
accession[is.na(accession) | !nzchar(accession)] <- NA_character_
symbol <- unname(gb2sym(accession))
coverage <- mean(!is.na(symbol) & nzchar(symbol))
cat(sprintf(
  "GSE56699: 58 pretreatment arrays -> 52 patients; 57 labeled arrays -> 51 labeled patients\n"
))
cat(sprintf(
  "GSE56699 patient groups: good=%d poor=%d; annotation coverage=%.1f%%\n",
  sum(grp == "good"), sum(grp == "poor"), 100 * coverage
))
stopifnot(coverage > 0.90)

y_sym <- symbolize(patient_y, symbol)
patient_features <- build_cohort_features(
  y_sym, grp, "GSE56699", ifn_sets, mcp_sets,
  write = TRUE, verbose = TRUE
)

# Sensitivity check for the repeated-biopsy collapse rule: score every array
# first and then average scores within patient.  The registered feature table
# above keeps the prespecified expression-average-then-score approach; this
# second route is an audit only and must not silently replace it.
array_y_sym <- symbolize(y, symbol)
array_grp <- array_response[keep]
names(array_grp) <- selected_assays
array_features <- build_cohort_features(
  array_y_sym, array_grp, "GSE56699_array_level", ifn_sets, mcp_sets,
  write = FALSE, verbose = FALSE
)
array_features$patient_code <- selected_patient[
  match(rownames(array_features), selected_assays)
]
core_features <- c("IFN", "NK", "Cyto")
score_average <- aggregate(
  array_features[core_features],
  list(patient_code = array_features$patient_code),
  mean
)
score_average <- score_average[
  match(rownames(patient_features), score_average$patient_code), , drop = FALSE
]
stopifnot(identical(score_average$patient_code, rownames(patient_features)))
score_average$response <- unname(grp[score_average$patient_code])

collapse_audit <- do.call(rbind, lapply(core_features, function(feature) {
  primary_g <- .hedges_g_soft(patient_features[[feature]], patient_features$response)
  alternate_g <- .hedges_g_soft(score_average[[feature]], score_average$response)
  data.frame(
    feature = feature,
    method = c(
      "mean_log2_expression_then_score",
      "score_each_array_then_mean"
    ),
    hedges_g = c(primary_g, alternate_g),
    correlation_with_primary = c(
      1,
      cor(patient_features[[feature]], score_average[[feature]], method = "spearman")
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  collapse_audit,
  "GSE56699_replicate_collapse_sensitivity.csv",
  row.names = FALSE
)
print(collapse_audit, row.names = FALSE, digits = 3)

cat("===== GSE56699 patient-level processing complete =====\n")
