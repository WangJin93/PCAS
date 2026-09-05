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

#' Resolve the on-disk cache directory
#'
#' Resolution order: explicit \code{cache_dir} argument > option
#' \code{PCAS.cache.dir} > the per-user cache directory
#' (\code{tools::R_user_dir("PCAS", "cache")}, i.e. \code{~/.cache/PCAS} on
#' Linux). Set the option to \code{NA} or \code{FALSE}, or pass
#' \code{use_cache = FALSE}, to disable caching entirely.
#' @noRd
.pcas_resolve_cache_dir <- function(use_cache = TRUE, cache_dir = NULL) {
  if (!isTRUE(use_cache)) return(NULL)
  if (is.character(cache_dir) && length(cache_dir) == 1L && nzchar(cache_dir)) {
    base <- cache_dir
  } else {
    cfg <- getOption("PCAS.cache.dir")
    if (identical(cfg, FALSE) || identical(cfg, NA)) return(NULL)
    base <- if (is.character(cfg) && length(cfg) == 1L && nzchar(cfg)) {
      cfg
    } else {
      tools::R_user_dir("PCAS", which = "cache")
    }
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

#' Compute the cache file name for one dataset/query (GCAS-style:
#' <dataset>_<md5 of ids>.RData under the <what> sub-directory)
#' @noRd
.pcas_cache_file <- function(cache_dir, dataset, ids, what = "data_temp") {
  h <- digest::digest(sort(unique(ids)), algo = "md5")
  file.path(cache_dir, what, paste0(dataset, "_", h, ".RData"))
}

#' Ensure a directory exists for a cache file path; returns path or NULL
#' @noRd
.pcas_ensure_cache_subdir <- function(file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(dirname(file))) file else NULL
}

#' Load a .RData cache file written by .pcas_save_cache() (object name
#' `cached_data`, as in the GCAS package); returns NULL on any failure.
#' @noRd
.pcas_load_cache <- function(file) {
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({
    load(file, envir = e)
    exists("cached_data", envir = e, inherits = FALSE)
  }, error = function(e2) FALSE)
  if (ok) e$cached_data else NULL
}

#' Save an object to a .RData cache file under the name `cached_data`;
#' best effort - never throws.
#' @noRd
.pcas_save_cache <- function(obj, file) {
  cached_data <- obj
  tryCatch({
    save(cached_data, file = file)
    TRUE
  }, error = function(e) {
    .pcas_note("Could not write cache file ", file,
               " (", conditionMessage(e), "); continuing without it.")
    FALSE
  })
}

#' Load a single-object .rds cache file (best effort)
#' @noRd
.pcas_load_rds <- function(file) {
  tryCatch(readRDS(file), error = function(e) NULL)
}

#' Save a single object as .rds (best effort)
#' @noRd
.pcas_save_rds <- function(obj, file) {
  tryCatch({
    saveRDS(obj, file)
    TRUE
  }, error = function(e) {
    .pcas_note("Could not write cache file ", file,
               " (", conditionMessage(e), "); continuing without it.")
    FALSE
  })
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
