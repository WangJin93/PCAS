#' @title Volcano plot for DEGs
#' @description
#' Volcano plot of differential expression results (from
#' \code{\link{get_DEGs_result}()}): x = log2 fold change, y = -log10
#' (adjusted p-value). Genes can be labelled by top up/down regulators or by a
#' user-supplied list.
#' @param df DEGs result obtained from \code{\link{get_DEGs_result}()}, i.e. a
#'   data.frame with a gene column and columns \code{logFC}, \code{P.Value}
#'   and/or \code{adj.P.Val}.
#' @param p.cut Adjusted p-value threshold (0 < p.cut < 1), default 0.05.
#' @param logFC.cut |log2 fold change| threshold, default 1.
#' @param show.top Label the 5 most down- and the 5 most up-regulated genes
#'   (needs >= 10 rows; otherwise labelling is skipped with a message).
#' @param show.labels Character vector of gene symbols to label; takes
#'   precedence over \code{show.top} when both are given.
#' @param label.size Size of the gene labels (default 5).
#' @return A ggplot object. Missing/insufficient columns and invalid thresholds
#'   are reported with messages instead of cryptic errors.
#' @examples
#' \dontrun{
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "limma")
#' viz_DEGs_volcano(results)
#' }
#' @export
viz_DEGs_volcano <- function(df,
                             p.cut = 0.05,
                             logFC.cut = 1,
                             show.top = FALSE,
                             show.labels = NULL,
                             label.size = 5) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    .pcas_note("viz_DEGs_volcano(): 'df' is empty; returning NULL.")
    return(NULL)
  }
  need <- c("logFC", "P.Value", "adj.P.Val")
  miss <- setdiff(need, colnames(df))
  if ("logFC" %in% miss) {
    stop("viz_DEGs_volcano(): input is missing the required 'logFC' column. ",
         call. = FALSE)
  }
  if ("adj.P.Val" %in% miss) {
    .pcas_note("No 'adj.P.Val' column found; falling back to 'P.Value' for ",
               "the y axis.")
  }

  p.cut <- suppressWarnings(as.numeric(p.cut))
  logFC.cut <- suppressWarnings(as.numeric(logFC.cut))
  if (is.na(p.cut) || p.cut <= 0 || p.cut >= 1) {
    warning("viz_DEGs_volcano(): invalid p.cut (", p.cut,
            "); using 0.05.", call. = FALSE)
    p.cut <- 0.05
  }
  if (is.na(logFC.cut) || logFC.cut < 0) {
    warning("viz_DEGs_volcano(): invalid logFC.cut (", logFC.cut,
            "); using 1.", call. = FALSE)
    logFC.cut <- 1
  }

  # label column: prefer "Symbol", else the first column
  lab_col <- if ("Symbol" %in% colnames(df)) "Symbol" else colnames(df)[1L]

  df$logFC <- .pcas_as_numeric(df$logFC)
  p_col <- if ("adj.P.Val" %in% colnames(df)) "adj.P.Val" else "P.Value"
  df[[p_col]] <- .pcas_as_numeric(df[[p_col]])

  na_rows <- is.na(df$logFC) | is.na(df[[p_col]])
  if (any(na_rows)) {
    .pcas_note("Removed ", sum(na_rows), " row(s) with missing logFC/p-value ",
               "before plotting.")
    df <- df[!na_rows, , drop = FALSE]
  }
  if (!nrow(df)) {
    .pcas_note("viz_DEGs_volcano(): no usable row left.")
    return(NULL)
  }

  yvals <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$neg_log10p <- pmin(yvals, 300)   # avoid Inf for p = 0

  df$change <- ifelse(df[[p_col]] < p.cut,
                      ifelse(df$logFC < -logFC.cut, "Down",
                             ifelse(df$logFC > logFC.cut, "Up", "No")),
                      "No")
  df$change <- factor(df$change, levels = c("Down", "No", "Up"))
  df$label <- ""

  if (!is.null(show.labels)) {
    show.labels <- unique(as.character(show.labels))
    found <- intersect(show.labels, as.character(df[[lab_col]]))
    missing_lab <- setdiff(show.labels, as.character(df[[lab_col]]))
    if (length(missing_lab)) {
      .pcas_note("Symbol(s) requested for labelling not present in the data: ",
                 paste(missing_lab, collapse = ", "))
    }
    df$label[as.character(df[[lab_col]]) %in% found] <-
      as.character(df[[lab_col]])[as.character(df[[lab_col]]) %in% found]
    if (show.top) {
      .pcas_note("Both show.top and show.labels are set; show.labels takes ",
                 "precedence.")
    }
  } else if (show.top) {
    if (nrow(df) < 10L) {
      .pcas_note("show.top = TRUE needs >= 10 rows (got ", nrow(df),
                 "); top-gene labelling skipped.")
    } else {
      o <- order(df$logFC)
      top_ids <- c(o[seq_len(5)], o[(nrow(df) - 4):nrow(df)])
      top_ids <- unique(top_ids)
      df$label[top_ids] <- as.character(df[[lab_col]])[top_ids]
    }
  }

  p <- ggplot2::ggplot(df) +
    ggplot2::aes(x = logFC, y = neg_log10p, colour = change) +
    ggplot2::geom_point(size = 2, shape = 16, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(Down = "blue", No = "grey20",
                                           Up = "red")) +
    ggplot2::geom_hline(yintercept = -log10(p.cut),
                        linetype = "dashed", colour = "grey30", size = 0.8) +
    ggplot2::geom_vline(xintercept = c(logFC.cut, -logFC.cut),
                        linetype = "dashed", colour = "grey30", size = 0.8) +
    ggplot2::theme_light(base_size = 20) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(
      x = bquote(''*Log[2]*' (Fold change)'),
      y = bquote(''*-Log[10]*' (P.adj value)'))

  if (any(nzchar(df$label))) {
    p <- p + ggrepel::geom_text_repel(
      data = df[nzchar(df$label), , drop = FALSE],
      ggplot2::aes(label = label),
      size = label.size, colour = "black",
      nudge_x = 0.2, nudge_y = 0.2,
      max.overlaps = 100000,
      box.padding = ggplot2::unit(0.9, "lines"),
      point.padding = ggplot2::unit(0.8, "lines"),
      show.legend = FALSE)
  }
  p
}
