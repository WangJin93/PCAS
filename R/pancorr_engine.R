# ============================================================================
# Internal engine shared by cor_pancancer_genelist / cor_pancancer_drug /
# cor_pancancer_TIL. NOT exported.
#
# Given one data.frame that already contains per-sample keys (ID, type,
# dataset), a target column and a set of feature columns, compute per-dataset
# correlation matrices r, p and a per-cell sample-size matrix n. Missing
# results stay NA (they are NOT silently removed with na.omit) and every
# dataset that was skipped for lack of samples is reported in the summary.
# ============================================================================

#' @noRd
.pancorr_core <- function(merged,
                          target_col,
                          feature_cols,
                          method = "pearson",
                          min_n = 4) {

  .pcas_require_cols(merged, c("ID", "type", "dataset", target_col),
                     "correlation analysis")

  method <- match.arg(method, c("pearson", "spearman"))

  # --- which requested features actually exist in the merged table ----------
  present_feats <- intersect(feature_cols, colnames(merged))
  absent_feats  <- setdiff(feature_cols, colnames(merged))
  if (length(absent_feats)) {
    .pcas_note("Feature(s) without data in the merged table and skipped: ",
               paste(utils::head(absent_feats, 10), collapse = ", "),
               if (length(absent_feats) > 10) paste0(" (+", length(absent_feats) - 10, " more)") else "")
  }
  if (!length(present_feats)) {
    .pcas_note("None of the requested features have data; no correlation can ",
               "be computed.")
    return(NULL)
  }

  # --- numeric coercion of target and features ------------------------------
  merged[[target_col]] <- .pcas_as_numeric(merged[[target_col]])
  for (cl in present_feats) merged[[cl]] <- .pcas_as_numeric(merged[[cl]])

  # --- per-dataset split (grouped by the stripped cohort label) -------------
  sss <- split(merged, merged$dataset)
  ds_labels <- names(sss)
  if (anyDuplicated(ds_labels)) {
    .pcas_note("Several input datasets map to the same cohort label; using ",
               "their full names in the output.")
    ds_labels <- vapply(sss, function(d) as.character(d$dataset[1L]),
                        character(1))
  }

  nr <- length(ds_labels)
  nc <- length(present_feats)
  r <- p <- n <- matrix(NA_real_, nrow = nr, ncol = nc,
                        dimnames = list(ds_labels, present_feats))

  skip_info <- data.frame(dataset = ds_labels, n_samples = NA_integer_,
                          reason = "", stringsAsFactors = FALSE)

  for (i in seq_along(sss)) {
    sub <- sss[[i]]
    skip_info$n_samples[i] <- nrow(sub)
    if (nrow(sub) < min_n) {
      skip_info$reason[i] <- sprintf(
        "only %d sample(s) (< min_n = %d)", nrow(sub), min_n)
      next
    }
    x <- sub[[target_col]]
    if (all(is.na(x))) {
      skip_info$reason[i] <- "target column has no value in this dataset"
      next
    }
    Y <- as.data.frame(lapply(sub[, present_feats, drop = FALSE],
                              as.numeric), check.names = FALSE)

    # per-feature complete-pair sample sizes
    for (j in seq_len(nc)) {
      n[i, j] <- sum(stats::complete.cases(x, Y[[j]]))
    }

    ct <- tryCatch(
      psych::corr.test(x = x, y = Y, method = method, ci = FALSE),
      error = function(e) e
    )
    if (inherits(ct, "error")) {
      skip_info$reason[i] <- paste("corr.test failed:",
                                   conditionMessage(ct))
      next
    }
    if (is.null(ct$r) || is.null(ct$p)) {
      skip_info$reason[i] <- "corr.test returned no estimate"
      next
    }
    rr <- as.numeric(ct$r); pp <- as.numeric(ct$p)
    r[i, ] <- rr
    p[i, ] <- pp
  }

  skipped <- skip_info[skip_info$reason != "" & !is.na(skip_info$reason), ,
                       drop = FALSE]
  if (nrow(skipped)) {
    .pcas_note("Dataset(s) excluded from correlation matrices: ",
               paste(sprintf("%s (%s)", skipped$dataset, skipped$reason),
                     collapse = "; "), ".")
  }
  if (sum(!is.na(r)) == 0L) {
    .pcas_note("No correlation could be computed in any dataset ",
               "(all matrices are NA).")
  }

  # --- keep a clean per-dataset table for the scatter-plot drill-down --------
  keep <- unique(c("ID", "type", "dataset", target_col, present_feats))
  keep <- intersect(keep, colnames(merged))
  sss <- lapply(sss, function(d) d[, keep, drop = FALSE])

  list(r = r, p = p, n = n, sss = sss,
       summary = list(
         target       = target_col,
         features     = present_feats,
         absent_features = absent_feats,
         min_n        = min_n,
         skipped      = skipped,
         datasets     = ds_labels,
         method       = method
       ))
}
