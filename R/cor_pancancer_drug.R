#' @title Correlation analysis of drug sensitivity
#' @description
#' Calculate the correlation between target gene expression and the sensitivity
#' of anti-tumour drugs (drug-response tables of the CPTAC pharmacoproteomics
#' data) in multiple datasets, dataset by dataset.
#' @param df The expression data of the target gene in multiple datasets,
#'   obtained by \code{\link{get_expr_data}()}. Only the first expression
#'   column is used as the target (a message is printed if more are supplied).
#' @param cor_method Correlation method, \code{"pearson"} or \code{"spearman"}.
#' @param Target.pathway The signalling pathways of the drug targets; use
#'   \code{unique(drug_info$Target.pathway)} for all valid values. Default
#'   \code{"Cell cycle"}.
#' @param min_n Minimum number of samples per dataset to run the correlation
#'   (default 4).
#' @return A list with elements \code{r}, \code{p}, \code{n} (dataset x drug
#'   matrices of correlation, p-value and pairwise-complete sample size),
#'   \code{sss} (per-dataset data.frames for scatter plots) and \code{summary}.
#'   \code{NULL} with an explanatory message when the pathway/features are
#'   invalid or nothing can be merged. Missing cells stay NA (never removed).
#' @details
#'   Only tumour samples that also have drug-sensitivity measurements can be
#'   used; the number of matched/mismatched samples per dataset is reported in
#'   \code{summary} and printed.
#' @examples
#' \dontrun{
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein",
#'              "LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein",
#'              "UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' result <- cor_pancancer_drug(df, Target.pathway = "Cell cycle")
#' }
#' @export
cor_pancancer_drug <- function(df,
                               cor_method = "pearson",
                               Target.pathway = "Cell cycle",
                               min_n = 4) {
  if (is.null(df)) {
    .pcas_note("cor_pancancer_drug(): 'df' is NULL; nothing to correlate.")
    return(NULL)
  }
  .pcas_require_cols(df, c("ID", "type", "dataset"), "cor_pancancer_drug")
  meta   <- c("ID", "type", "dataset")
  target <- setdiff(colnames(df), meta)[1L]
  if (is.na(target)) {
    .pcas_note("cor_pancancer_drug(): 'df' has no expression column.")
    return(NULL)
  }
  expr_df <- setdiff(colnames(df), meta)
  if (length(expr_df) > 1L) {
    .pcas_note("'df' contains ", length(expr_df),
               " expression columns; using only the first one ('", target,
               "') as the correlation target.")
  }
  df$dataset <- .pcas_strip_dataset(df$dataset)

  # ---- validate the pathway argument with suggestions -----------------------
  if (!exists("drug_info") || !is.data.frame(drug_info)) {
    stop("cor_pancancer_drug(): the drug_info object is not available.",
         call. = FALSE)
  }
  valid_pathways <- unique(drug_info$Target.pathway)
  bad <- setdiff(Target.pathway, valid_pathways)
  if (length(bad)) {
    suggestion <- .pcas_suggest(bad, valid_pathways)
    warning("cor_pancancer_drug(): unknown Target.pathway value(s): ",
            paste(bad, collapse = ", "),
            paste0(" (did you mean '",
                   paste(suggestion[!is.na(suggestion)], collapse = "', '"), "'?) "),
            "Valid choices: ", paste(valid_pathways, collapse = ", "), ".",
            call. = FALSE)
    Target.pathway <- setdiff(Target.pathway, bad)
    if (!length(Target.pathway)) return(NULL)
  }

  # ---- drugs of the selected pathway that are present in the data -----------
  sub_info <- drug_info[drug_info$Target.pathway %in% Target.pathway,
                        c("ID", "Name")]
  sig  <- unique(paste(sub_info$Name, sub_info$ID, sep = "_"))
  miss <- setdiff(sig, colnames(drug_CPTAC))
  if (length(miss)) {
    .pcas_note("Drug(s) annotated in drug_info but absent from the drug ",
               "sensitivity matrix and skipped: ",
               paste(utils::head(miss, 10), collapse = ", "),
               if (length(miss) > 10) paste0(" (+", length(miss) - 10,
                                             " more)") else "")
  }
  sig <- intersect(sig, colnames(drug_CPTAC))
  if (!length(sig)) {
    .pcas_note("No drug of the selected pathway has sensitivity data; ",
               "nothing to correlate.")
    return(NULL)
  }

  # ---- merge expression with drug sensitivity -------------------------------
  drug_df <- drug_CPTAC[, sig, drop = FALSE]
  drug_df <- tibble::rownames_to_column(as.data.frame(drug_df), var = "ID")

  pre_counts <- table(df$dataset)
  merged <- merge(df, drug_df, by = "ID", all = FALSE)
  post_counts <- table(merged$dataset)
  for (ds in names(pre_counts)) {
    before <- pre_counts[[ds]]
    after  <- if (ds %in% names(post_counts)) post_counts[[ds]] else 0L
    if (after < before) {
      .pcas_note("Dataset ", ds, ": ", before - after, " of ", before,
                 " sample row(s) had no drug-sensitivity measurement and were ",
                 "left out of the correlation.")
    }
  }
  if (!nrow(merged)) {
    .pcas_note("cor_pancancer_drug(): no sample matched the drug-sensitivity ",
               "table.")
    return(NULL)
  }

  core <- .pancorr_core(merged, target, sig, method = cor_method,
                        min_n = min_n)
  if (is.null(core)) return(NULL)
  core$summary$Target.pathway <- Target.pathway
  list(r = core$r, p = core$p, n = core$n, sss = core$sss,
       summary = core$summary)
}
