#' @title Perform pan-cancer correlation analysis against a gene set
#' @description
#' Correlate the expression of one target gene with a list of genes/features
#' across several CPTAC datasets, dataset by dataset. Returns correlation and
#' p-value matrices as well as a sample-size matrix and per-dataset data for
#' scatter plots.
#' @param df The expression data of the target gene in multiple datasets,
#'   obtained by \code{\link{get_expr_data}()}. Only the first expression
#'   column is used as the target (a message is printed if more are supplied).
#' @param geneset_data The expression data of a genelist in multiple datasets,
#'   obtained by \code{\link{get_expr_data}()}.
#' @param sample_type Sample type used for correlation analysis, e.g.
#'   \code{c("Tumor", "Normal")}.
#' @param cor_method Correlation method, \code{"pearson"} or \code{"spearman"}.
#' @param min_n Minimum number of samples per dataset to run the correlation
#'   (datasets below this are skipped and reported, not silently dropped).
#' @return A list with elements \code{r}, \code{p}, \code{n} (dataset x feature
#'   matrices of correlation, p-value and pairwise-complete sample size),
#'   \code{sss} (per-dataset data.frames for scatter plots; columns
#'   ID/type/dataset/target/features) and \code{summary}. \code{NULL} with an
#'   explanatory message when nothing can be computed. Missing cells stay NA.
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7",
#'               "HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("CCRCC_CPTAC_protein","GBM_CPTAC_protein","HNSCC_CPTAC_protein",
#'              "LSCC_CPTAC_protein","LUAD_CPTAC_protein","PDAC_CPTAC_protein",
#'              "UCEC_CPTAC2_protein","UCEC_CPTAC1_protein")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist, datasets = dataset)
#' result <- cor_pancancer_genelist(df, geneset_data, sample_type = "Tumor")
#' }
#' @export
cor_pancancer_genelist <- function(df,
                                   geneset_data,
                                   sample_type = c("Tumor", "Normal"),
                                   cor_method = "pearson",
                                   min_n = 4) {
  if (is.null(df) || is.null(geneset_data)) {
    .pcas_note("cor_pancancer_genelist(): 'df'/'geneset_data' is NULL; ",
               "nothing to correlate.")
    return(NULL)
  }
  .pcas_require_cols(df, c("ID", "type", "dataset"),
                     "cor_pancancer_genelist")
  .pcas_require_cols(geneset_data, c("ID", "type", "dataset"),
                     "cor_pancancer_genelist")

  meta   <- c("ID", "type", "dataset")
  target <- setdiff(colnames(df), meta)[1L]
  if (is.na(target)) {
    .pcas_note("cor_pancancer_genelist(): 'df' has no expression column.")
    return(NULL)
  }
  expr_df <- setdiff(colnames(df), meta)
  if (length(expr_df) > 1L) {
    .pcas_note("'df' contains ", length(expr_df),
               " expression columns; using only the first one ('", target,
               "') as the correlation target.")
  }
  feature_cols <- setdiff(colnames(geneset_data), meta)
  if (!length(feature_cols)) {
    .pcas_note("cor_pancancer_genelist(): 'geneset_data' has no expression ",
               "column.")
    return(NULL)
  }
  sample_type <- intersect(sample_type, c("Tumor", "Normal", "Other"))
  if (!length(sample_type)) {
    warning("cor_pancancer_genelist(): 'sample_type' must be a subset of ",
            "c('Tumor', 'Normal', 'Other').", call. = FALSE)
    return(NULL)
  }

  # strip the _protein/_mRNA/_Phospho suffix so protein and mRNA tables of one
  # cohort are merged with each other
  df$dataset          <- .pcas_strip_dataset(df$dataset)
  geneset_data$dataset <- .pcas_strip_dataset(geneset_data$dataset)

  d1 <- df[df$type %in% sample_type, c(meta, target), drop = FALSE]
  d2 <- geneset_data[geneset_data$type %in% sample_type, , drop = FALSE]

  merged <- merge(d1, d2, by = meta, all = FALSE)
  n_unmatched <- nrow(d1) - nrow(merged)
  if (n_unmatched > 0L) {
    .pcas_note(n_unmatched, " target-gene sample row(s) (of ", nrow(d1),
               ") could not be matched with the gene-set data and were left ",
               "out (per-dataset sample sizes are in the 'n' matrix).")
  }
  if (!nrow(merged)) {
    .pcas_note("cor_pancancer_genelist(): no sample row could be matched ",
               "between 'df' and 'geneset_data'.")
    return(NULL)
  }

  core <- .pancorr_core(merged, target, feature_cols,
                        method = cor_method, min_n = min_n)
  if (is.null(core)) return(NULL)
  core$summary$sample_type <- sample_type
  list(r = core$r, p = core$p, n = core$n, sss = core$sss,
       summary = core$summary)
}
