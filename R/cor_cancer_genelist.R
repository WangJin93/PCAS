#' @title Perform single-cohort correlation analysis
#' @description
#' Correlate one identifier (gene / phosphosite) of one dataset with one or
#' several identifiers of a second dataset (protein, mRNA or phospho table of
#' the same cohort), on the samples shared by both tables.
#' @param dataset1 Dataset name of the first (target) identifier. Use
#'   \code{dataset_info$Abbre} to list all datasets.
#' @param id1 Single gene symbol or phosphosite id of \code{dataset1}.
#' @param dataset2 Dataset name of the second (feature) identifiers.
#' @param id2 One or several gene symbols / phosphosite ids of
#'   \code{dataset2}.
#' @param sample_type Sample types used, default \code{c("Tumor", "Normal")}.
#' @param cor_method Correlation method, \code{"pearson"} or \code{"spearman"}.
#' @param min_n Minimum number of shared samples to compute a correlation.
#' @return A list with elements
#'   \describe{
#'     \item{\code{cor_result}}{data.frame with columns \code{Symbol},
#'       \code{Correlation}, \code{P value} and \code{n} — one row per entry of
#'       \code{id2} (same order, missing features are reported and stay NA).}
#'     \item{\code{cor_data}}{the merged long-format data used for the test
#'       (columns ID/type/dataset/target/features), ready for scatter plots.}
#'     \item{\code{n}}{vector of pairwise-complete sample sizes per feature.}
#'     \item{\code{summary}}{matching and missingness report.}
#'   }
#'   Returns \code{NULL} (with messages) when no data could be fetched or
#'   matched.
#' @examples
#' \dontrun{
#' results <- cor_cancer_genelist(dataset1 = "LUAD_CPTAC_protein",
#'                                id1 = "STAT3",
#'                                dataset2 = "LUAD_CPTAC_mRNA",
#'                                id2 = c("TNS1", "TP53"),
#'                                sample_type = c("Tumor", "Normal"),
#'                                cor_method = "pearson")
#' }
#' @export
cor_cancer_genelist <- function(dataset1 = "LUAD_CPTAC_protein",
                                id1 = "STAT3",
                                dataset2 = "LUAD_CPTAC_mRNA",
                                id2 = c("TNS1", "TP53"),
                                sample_type = c("Tumor", "Normal"),
                                cor_method = "pearson",
                                min_n = 4) {
  meta <- c("ID", "type", "dataset")

  # ---------------------------------------------------------------------------
  # fetch both tables (NULL-check BEFORE any transformation)
  # ---------------------------------------------------------------------------
  data1 <- get_expr_data(dataset1, id1)
  if (is.null(data1)) {
    .pcas_note("cor_cancer_genelist(): no expression data for dataset1 '",
               dataset1, "' / id1 '", paste(id1, collapse = ","),
               "'; returning NULL.")
    return(NULL)
  }
  data2 <- get_expr_data(dataset2, id2)
  if (is.null(data2)) {
    .pcas_note("cor_cancer_genelist(): no expression data for dataset2 '",
               dataset2, "' / id2 '", paste(id2, collapse = ","),
               "'; returning NULL.")
    return(NULL)
  }
  .pcas_require_cols(data1, meta, "cor_cancer_genelist")
  .pcas_require_cols(data2, meta, "cor_cancer_genelist")

  data1$dataset <- .pcas_strip_dataset(data1$dataset)
  data2$dataset <- .pcas_strip_dataset(data2$dataset)

  sample_type <- intersect(sample_type, c("Tumor", "Normal", "Other"))
  if (!length(sample_type)) {
    warning("cor_cancer_genelist(): 'sample_type' must be a subset of ",
            "c('Tumor', 'Normal', 'Other').", call. = FALSE)
    return(NULL)
  }

  # target = first expression column of data1
  expr1 <- setdiff(colnames(data1), meta)
  target <- expr1[1L]
  if (length(expr1) > 1L) {
    .pcas_note("data1 contains ", length(expr1),
               " expression columns; using '", target, "' as the target.")
  }
  # features = expression columns of data2, in the order requested in id2
  id2 <- unique(id2)
  feat_cols <- intersect(id2, setdiff(colnames(data2), meta))
  missing_id2 <- setdiff(id2, feat_cols)
  if (length(missing_id2)) {
    .pcas_note("Identifier(s) of id2 with no measurement in dataset2: ",
               paste(missing_id2, collapse = ", "),
               " (rows are kept in the result but stay NA).")
  }

  d1 <- data1[data1$type %in% sample_type,
              unique(c(meta, target)), drop = FALSE]
  d2 <- data2[data2$type %in% sample_type, , drop = FALSE]

  merged <- merge(d1, d2, by = meta, all = FALSE)
  n_unmatched <- nrow(d1) - nrow(merged)
  if (n_unmatched > 0L) {
    .pcas_note(n_unmatched, " sample row(s) of dataset1 (of ", nrow(d1),
               ") are not measured in dataset2 and were left out of the ",
               "correlation.")
  }
  if (!nrow(merged)) {
    .pcas_note("cor_cancer_genelist(): no shared sample between dataset1 and ",
               "dataset2 for the requested sample types.")
    return(NULL)
  }
  # keep missing features as all-NA columns so result rows align with id2
  for (g in setdiff(id2, feat_cols)) merged[[g]] <- NA_real_

  x <- .pcas_as_numeric(merged[[target]])
  ycols <- id2
  Y <- as.data.frame(lapply(merged[, ycols, drop = FALSE],
                            .pcas_as_numeric), check.names = FALSE)

  nvec <- vapply(Y, function(y) sum(stats::complete.cases(x, y)), integer(1))

  ct <- tryCatch(psych::corr.test(x = x, y = Y, method = cor_method,
                                  ci = FALSE),
                 error = function(e) e)
  if (inherits(ct, "error") || is.null(ct$r)) {
    warning("cor_cancer_genelist(): correlation test failed for dataset1 '",
            dataset1, "' vs dataset2 '", dataset2, "': ",
            if (inherits(ct, "error")) conditionMessage(ct) else "no estimate",
            call. = FALSE)
    return(NULL)
  }
  rr <- as.numeric(ct$r); pp <- as.numeric(ct$p)

  cor_result <- data.frame(Symbol = ycols,
                           Correlation = rr,
                           "P value" = pp,
                           n = nvec,
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
  rownames(cor_result) <- NULL

  cor_data <- merged[, unique(c(meta, target, ycols)), drop = FALSE]

  summary <- list(
    dataset1 = dataset1, id1 = target,
    dataset2 = dataset2, id2 = id2,
    sample_type = sample_type,
    cor_method = cor_method,
    n_shared_samples = nrow(merged),
    n_unmatched_samples = n_unmatched,
    missing_id2 = missing_id2,
    n = nvec
  )
  if (nrow(merged) < min_n) {
    .pcas_note("Only ", nrow(merged), " shared sample(s) (< min_n = ", min_n,
               "); correlations may be unreliable.")
  }

  list(cor_result = cor_result, cor_data = cor_data, n = nvec,
       summary = summary)
}
