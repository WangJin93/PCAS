#' @title Heat map of correlation results
#' @description
#' Present correlation results (from \code{\link{cor_pancancer_genelist}()},
#' \code{\link{cor_pancancer_TIL}()} or \code{\link{cor_pancancer_drug}()}) as
#' a heat map: colour = correlation coefficient, text = significance stars.
#' Cells without data (NA, e.g. too few samples) keep their own grey colour and
#' are counted in a message — they are never turned into a coefficient of zero.
#' @param r The correlation coefficient matrix of the results.
#' @param p The p-value matrix of the results (same dimensions/rownames).
#' @param limits Colour scale limits, default \code{c(-1, 1)}.
#' @return A ggplot object. A hierarchical cluster tree is added on the left
#'   when the \pkg{ggtree}/\pkg{aplot} packages are installed and the matrix
#'   has no missing value (otherwise a message explains the fallback).
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7",
#'               "HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <- c("CCRCC_CPTAC_mRNA","GBM_CPTAC_mRNA","HNSCC_CPTAC_mRNA",
#'              "LSCC_CPTAC_mRNA","LUAD_CPTAC_mRNA","PDAC_CPTAC_mRNA",
#'              "UCEC_CPTAC2_mRNA","UCEC_CPTAC1_mRNA")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist, datasets = dataset)
#' result <- cor_pancancer_genelist(df, geneset_data, sample_type = "Tumor")
#' viz_cor_heatmap(result$r, result$p)
#' }
#' @export
viz_cor_heatmap <- function(r, p, limits = c(-1, 1)) {
  if (is.null(r) || is.null(p)) {
    .pcas_note("viz_cor_heatmap(): 'r'/'p' is NULL.")
    return(NULL)
  }
  r <- as.matrix(r); p <- as.matrix(p)
  if (!identical(dim(r), dim(p))) {
    stop("viz_cor_heatmap(): 'r' and 'p' must have the same dimensions.",
         call. = FALSE)
  }
  if (!identical(rownames(r), rownames(p)) ||
      !identical(colnames(r), colnames(p))) {
    stop("viz_cor_heatmap(): 'r' and 'p' must share the same row/column names.",
         call. = FALSE)
  }
  if (!nrow(r) || !ncol(r)) {
    .pcas_note("viz_cor_heatmap(): empty matrix.")
    return(NULL)
  }

  n_na_r <- sum(is.na(r))
  if (n_na_r > 0L) {
    .pcas_note("The correlation matrix contains ", n_na_r,
               " missing cell(s) (e.g. datasets with too few samples); they ",
               "are drawn in grey and are NOT treated as a correlation of 0.")
  }

  # ---- optional row clustering (needs ggtree + aplot) ------------------------
  row_order <- seq_len(nrow(r))
  tree_added <- FALSE
  if (nrow(r) > 1L) {
    has_pkgs <- requireNamespace("ggtree", quietly = TRUE) &&
      requireNamespace("aplot", quietly = TRUE)
    if (!has_pkgs) {
      .pcas_note("Packages 'ggtree' and 'aplot' are not installed; plotting ",
                 "the heat map without the cluster tree.")
    } else if (any(!is.finite(r))) {
      .pcas_note("The correlation matrix contains missing cells; the cluster ",
                 "tree is skipped.")
    } else {
      hc <- tryCatch(stats::hclust(stats::dist(r)), error = function(e) NULL)
      if (is.null(hc)) {
        .pcas_note("Row clustering failed; plotting without the tree.")
      } else {
        row_order <- hc$order
        tree_added <- TRUE
      }
    }
  }
  r <- r[row_order, , drop = FALSE]

  melt_corr <- reshape2::melt(as.matrix(r), na.rm = FALSE)
  melt_p    <- reshape2::melt(as.matrix(p), na.rm = FALSE)
  names(melt_corr) <- c("row", "col", "corr")
  names(melt_p)    <- c("row", "col", "pval")
  melt_data <- merge(melt_corr, melt_p, by = c("row", "col"),
                     sort = FALSE, all = TRUE)

  melt_data$text <- ifelse(is.na(melt_data$pval), "",
                           ifelse(melt_data$pval < 0.001, "***",
                                  ifelse(melt_data$pval < 0.01, "**",
                                         ifelse(melt_data$pval < 0.05, "*",
                                                ""))))
  # factor levels fix the plotting order of rows/columns
  col_clean <- stringr::str_remove(colnames(r), "_XCELL")
  melt_data$row <- factor(melt_data$row, levels = rownames(r))
  melt_data$col <- factor(melt_data$col, levels = col_clean)

  heat <- ggplot2::ggplot(melt_data, ggplot2::aes(col, row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = corr), colour = "grey",
                       linewidth = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#5C5DAF", mid = "white", high = "#EA2E2D",
      midpoint = 0, limits = limits, oob = scales::squish,
      na.value = "grey85", name = "Correlation") +
    ggplot2::geom_text(ggplot2::aes(label = text), colour = "black",
                       size = 7, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                          size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(l = 30, unit = "pt")) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::labs(fill = paste0("  *   p < 0.05", "\n",
                                " **  p < 0.01", "\n",
                                "*** p < 0.001", "\n\n", "Correlation"))

  if (tree_added) {
    hc <- stats::hclust(stats::dist(r))          # recompute with the final order
    tr <- ggtree::ggtree(hc, layout = "rectangular", branch.length = "none")
    heat <- aplot::insert_left(heat, tr, width = 0.1)
  }
  heat
}
