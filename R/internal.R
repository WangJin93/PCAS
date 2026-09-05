# ============================================================================
# Internal helpers for PCAS (not exported)
# ============================================================================

#' Strip the molecular-type suffix from a dataset abbreviation
#'
#' `LUAD_CPTAC_protein`, `LUAD_CPTAC_mRNA` and `LUAD_CPTAC_Phospho` all map to
#' the cohort label `LUAD_CPTAC`.
#' @noRd
.pcas_strip_dataset <- function(x) {
  sub("_(protein|Phospho|mRNA)$", "", x)
}

#' Parse sample type and sample ID out of a CPTAC sample column name
#'
#' Sample columns look like `C3L-00094_Tumor`, `C3L-00094_Normal` or pooled
#' columns such as `Normal Only IR_Other`. Returns a data.frame with columns
#' `ID` and `type`.
#' @noRd
.pcas_parse_sample <- function(sample_names) {
  type <- sub("^.*_([^_]+)$", "\\1", sample_names)
  id   <- sub("_(Tumor|Normal|Other)$", "", sample_names)
  data.frame(ID = id, type = type, stringsAsFactors = FALSE)
}

#' Cache directory for PCAS queries
#'
#' Location is controlled with \code{options(PCAS.cache.dir = <path>)}.
#' Set it to \code{NA} or \code{FALSE} to disable the on-disk cache.
#' Defaults to a per-user cache directory (R >= 4.0).
#' @noRd
.pcas_cache_dir <- function(use_cache = TRUE) {
  if (!isTRUE(use_cache)) return(NULL)
  cfg <- getOption("PCAS.cache.dir")
  if (identical(cfg, FALSE) || identical(cfg, NA)) return(NULL)
  base <- if (is.character(cfg) && length(cfg) == 1L && nzchar(cfg)) {
    cfg
  } else {
    tools::R_user_dir("PCAS", which = "cache")
  }
  ok <- tryCatch({
    dir.create(base, recursive = TRUE, showWarnings = FALSE)
    dir.exists(base)
  }, error = function(e) FALSE)
  if (!ok) {
    warning("Cannot create PCAS cache directory '", base,
            "'; continuing without on-disk caching.", call. = FALSE)
    return(NULL)
  }
  base
}

#' Compute a stable cache file name for one dataset/query
#' @noRd
.pcas_cache_file <- function(cache_dir, dataset, ids, what = "expr") {
  h <- digest::digest(sort(unique(ids)), algo = "md5")
  month <- format(Sys.Date(), "%Y-%m")   # invalidates stale DB content monthly
  file.path(cache_dir, what, paste0(dataset, "_", h, "_", month, ".RData"))
}

#' Ensure a directory exists for a cache file path; returns path or NULL
#' @noRd
.pcas_ensure_cache_subdir <- function(file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(dirname(file))) file else NULL
}

#' Cast columns to numeric and report how many NAs are introduced
#' @noRd
.pcas_as_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(y <- as.numeric(x))
  new_na <- sum(is.na(y) & !is.na(x))
  if (new_na > 0L) {
    warning("Coercing a non-numeric column to numeric introduced ", new_na,
            " missing value(s).", call. = FALSE)
  }
  y
}

#' Unified, tidy status message
#' @noRd
.pcas_note <- function(...) message("PCAS: ", ...)

#' Validate that required columns exist in a data.frame
#' @noRd
.pcas_require_cols <- function(x, cols, fun) {
  miss <- setdiff(cols, colnames(x))
  if (length(miss)) {
    stop(fun, "(): input is missing required column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Recommend a similar valid value for a typo'd categorical input
#' @noRd
.pcas_suggest <- function(bad, valid, n = 3L) {
  if (!length(bad) || !length(valid)) return(character(0))
  out <- character(length(bad))
  for (i in seq_along(bad)) {
    d <- utils::adist(bad[i], valid)
    out[i] <- if (min(d) <= 2L) valid[which.min(d)] else NA_character_
  }
  out
}
