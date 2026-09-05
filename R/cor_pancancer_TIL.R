#' @title Correlation analysis of immune cell infiltration
#' @description
#' Calculate the correlation between target gene expression and immune cell
#' infiltration (tumour samples only) in multiple datasets, dataset by dataset.
#' @param df The expression data of the target gene in multiple datasets,
#'   obtained by \code{\link{get_expr_data}()}. Only the first expression
#'   column is used as the target (a message is printed if more are supplied).
#' @param cor_method Correlation method, \code{"pearson"} or \code{"spearman"}
#'   (default \code{"spearman"}).
#' @param TIL_type Algorithm(s) used to derive the infiltration scores, e.g.
#'   \code{"TIMER"}, \code{"MCPcounter"}, \code{"EPIC"} ... Use
#'   \code{unique(TIL_map$algorithm)} for all valid values.
#' @param min_n Minimum number of samples per dataset to run the correlation
#'   (default 4).
#' @return A list with elements \code{r}, \code{p}, \code{n} (dataset x cell
#'   type matrices of correlation, p-value and pairwise-complete sample size),
#'   \code{sss} (per-dataset data.frames for scatter plots) and \code{summary}.
#'   \code{NULL} with an explanatory message when the algorithm is invalid or
#'   nothing can be merged. Missing cells stay NA (never removed).
#' @examples
#' \dontrun{
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein",
#'              "LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein",
#'              "UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' result <- cor_pancancer_TIL(df, TIL_type = "TIMER")
#' }
#' @export
cor_pancancer_TIL <- function(df,
                              cor_method = "spearman",
                              TIL_type = "TIMER",
                              min_n = 4) {
  if (is.null(df)) {
    .pcas_note("cor_pancancer_TIL(): 'df' is NULL; nothing to correlate.")
    return(NULL)
  }
  .pcas_require_cols(df, c("ID", "type", "dataset"), "cor_pancancer_TIL")
  meta   <- c("ID", "type", "dataset")
  target <- setdiff(colnames(df), meta)[1L]
  if (is.na(target)) {
    .pcas_note("cor_pancancer_TIL(): 'df' has no expression column.")
    return(NULL)
  }
  expr_df <- setdiff(colnames(df), meta)
  if (length(expr_df) > 1L) {
    .pcas_note("'df' contains ", length(expr_df),
               " expression columns; using only the first one ('", target,
               "') as the correlation target.")
  }
  df$dataset <- .pcas_strip_dataset(df$dataset)

  # ---- validate the TIL algorithm argument with suggestions -----------------
  if (!exists("TIL_map") || !is.data.frame(TIL_map)) {
    stop("cor_pancancer_TIL(): the TIL_map object is not available.",
         call. = FALSE)
  }
  valid_algo <- unique(TIL_map$algorithm)
  bad <- setdiff(TIL_type, valid_algo)
  if (length(bad)) {
    suggestion <- .pcas_suggest(bad, valid_algo)
    warning("cor_pancancer_TIL(): unknown TIL_type value(s): ",
            paste(bad, collapse = ", "),
            paste0(" (did you mean '",
                   paste(suggestion[!is.na(suggestion)], collapse = "', '"), "'?) "),
            "Valid choices: ", paste(valid_algo, collapse = ", "), ".",
            call. = FALSE)
    TIL_type <- setdiff(TIL_type, bad)
    if (!length(TIL_type)) return(NULL)
  }

  sig <- TIL_map$cell_type[TIL_map$algorithm %in% TIL_type]
  miss <- setdiff(sig, colnames(TIL_CPTAC))
  if (length(miss)) {
    .pcas_note("Cell type(s) listed in TIL_map but absent from TIL_CPTAC and ",
               "skipped: ", paste(utils::head(miss, 10), collapse = ", "),
               if (length(miss) > 10) paste0(" (+", length(miss) - 10,
                                             " more)") else "")
  }
  sig <- intersect(sig, colnames(TIL_CPTAC))
  if (!length(sig)) {
    .pcas_note("No infiltration column available for the selected TIL ",
               "algorithm(s); nothing to correlate.")
    return(NULL)
  }

  # ---- tumour samples only, merged with infiltration scores -----------------
  n_normal <- sum(df$type != "Tumor")
  if (n_normal > 0L) {
    .pcas_note("Immune infiltration is scored on tumour samples only; ",
               n_normal, " non-tumour sample row(s) were not used.")
  }
  df_tumour <- df[df$type == "Tumor", , drop = FALSE]
  if (!nrow(df_tumour)) {
    .pcas_note("cor_pancancer_TIL(): no tumour sample in 'df'.")
    return(NULL)
  }
  til <- TIL_CPTAC[, c("ID", sig), drop = FALSE]
  til <- as.data.frame(til)
  merged <- merge(df_tumour, til, by = "ID", all = FALSE)
  if (!nrow(merged)) {
    .pcas_note("cor_pancancer_TIL(): no tumour sample matched the immune ",
               "infiltration table.")
    return(NULL)
  }

  core <- .pancorr_core(merged, target, sig, method = cor_method,
                        min_n = min_n)
  if (is.null(core)) return(NULL)
  core$summary$TIL_type <- TIL_type
  list(r = core$r, p = core$p, n = core$n, sss = core$sss,
       summary = core$summary)
}
