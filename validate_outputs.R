# Validate cohort-level feature tables before meta-analysis.
# This deliberately fails when a stale, duplicated or incorrectly labelled table
# is present, while reporting legacy extra columns as warnings for transparency.

if (!exists("LARC_CONFIG", inherits = FALSE)) {
  if (!file.exists("project_config.R")) {
    stop("project_config.R not found; run this script from the project directory")
  }
  source("project_config.R", encoding = "UTF-8")
}
cfg <- LARC_CONFIG
expected <- cfg$cohorts
required_columns <- cfg$required_feature_columns
numeric_columns <- cfg$meta_feature_columns

validate_one <- function(path, expected_row) {
  if (!file.exists(path)) stop("Missing cohort feature table: ", path)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dat) < 2L) stop(path, " must contain a sample identifier and feature columns")

  missing_columns <- setdiff(required_columns, colnames(dat))
  if (length(missing_columns)) {
    stop(path, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  extra_columns <- setdiff(colnames(dat)[-1L], required_columns)
  if (length(extra_columns)) {
    warning(path, " contains extra columns: ", paste(extra_columns, collapse = ", "),
            ". They are ignored by the meta-analysis; regenerate to remove legacy fields.")
  }

  sample_id <- trimws(as.character(dat[[1L]]))
  if (anyNA(sample_id) || any(!nzchar(sample_id))) {
    stop(path, " contains empty sample identifiers")
  }
  if (anyDuplicated(sample_id)) {
    dup <- unique(sample_id[duplicated(sample_id)])
    stop(path, " contains duplicated sample identifiers: ", paste(head(dup, 5L), collapse = ", "))
  }

  response <- as.character(dat$response)
  if (anyNA(response) || any(!response %in% c("good", "poor"))) {
    stop(path, " has invalid response labels; expected only good/poor")
  }
  for (column in numeric_columns) {
    values <- suppressWarnings(as.numeric(dat[[column]]))
    if (any(!is.finite(values))) {
      stop(path, " contains missing or non-finite meta-analysis values in ", column)
    }
  }

  observed <- c(
    n = nrow(dat),
    n_good = sum(response == "good"),
    n_poor = sum(response == "poor")
  )
  target <- as.integer(unlist(expected_row[c("n", "n_good", "n_poor")],
                              use.names = FALSE))
  if (!identical(unname(as.integer(observed)), unname(target))) {
    stop(
      path, " has stale or unexpected counts. Observed ",
      paste(names(observed), observed, sep = "=", collapse = ", "),
      "; expected ", paste(names(target), target, sep = "=", collapse = ", ")
    )
  }

  data.frame(
    cohort = expected_row$cohort,
    config_version = cfg$version,
    n = unname(observed["n"]),
    n_good = unname(observed["n_good"]),
    n_poor = unname(observed["n_poor"]),
    sample_id_unique = TRUE,
    extra_columns = if (length(extra_columns)) paste(extra_columns, collapse = ";") else "",
    check = "PASS",
    stringsAsFactors = FALSE
  )
}

checks <- lapply(seq_len(nrow(expected)), function(i) {
  cohort <- expected$cohort[i]
  validate_one(paste0("features_", cohort, ".csv"), expected[i, , drop = FALSE])
})
checks <- do.call(rbind, checks)
rownames(checks) <- NULL

stopifnot(
  sum(checks$n) == sum(cfg$cohorts$n),
  sum(checks$n[checks$cohort %in% cfg$analysis_sets$primary_ncrt9]) ==
    sum(cfg$cohorts$n[cfg$cohorts$treatment_scope == "nCRT"]),
  all(checks$check == "PASS")
)
write.csv(checks, "validation_summary.csv", row.names = FALSE)
print(checks, row.names = FALSE)
cat(paste0(
  "Validated: 9 nCRT cohorts (n=577), 1 radiotherapy-only sensitivity cohort ",
  "(n=51 independent patients).\n"
))
