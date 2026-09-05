#' @title Scatter plot of two variables with correlation annotation
#' @description
#' Scatter plot of two expression columns with linear smooth, rug, and a title
#' giving the sample size, correlation coefficient and p-value
#' (\code{\link[stats]{cor.test}}).
#' @param data A data.frame containing the two columns \code{a} and \code{b}.
#' @param a Column name of variable A.
#' @param b Column name of variable B.
#' @param method Correlation method, \code{"pearson"} or \code{"spearman"}.
#' @param x_lab X-axis label suffix.
#' @param y_lab Y-axis label suffix.
#' @return A ggplot object. Rows with missing values are counted and reported,
#'   and the correlation is computed on the complete observations only.
#' @examples
#' \dontrun{
#' df <- get_expr_data("LUAD_CPTAC_protein", c("TP53", "TNS1"))
#' viz_corplot(df, "TP53", "TNS1")
#' }
#' @export
viz_corplot <- function(data,
                        a, b,
                        method = "pearson",
                        x_lab = " expression",
                        y_lab = " expression") {
  if (is.null(data) || !is.data.frame(data)) {
    .pcas_note("viz_corplot(): 'data' is NULL/not a data.frame.")
    return(NULL)
  }
  method <- match.arg(method, c("pearson", "spearman"))
  miss <- setdiff(c(a, b), colnames(data))
  if (length(miss)) {
    stop("viz_corplot(): column(s) not found in 'data': ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  data <- data[, c(a, b), drop = FALSE]
  data[[a]] <- .pcas_as_numeric(data[[a]])
  data[[b]] <- .pcas_as_numeric(data[[b]])
  names(data) <- c("geneA", "geneB")

  n_na <- sum(!stats::complete.cases(data))
  if (n_na > 0L) {
    .pcas_note("Removed ", n_na, " row(s) with missing values (of ", nrow(data),
               ") before correlation/plotting.")
  }
  comp <- data[stats::complete.cases(data), , drop = FALSE]
  if (nrow(comp) < 3L) {
    .pcas_note("Fewer than 3 complete observations (", nrow(comp),
               "); correlation cannot be estimated reliably.")
  }

  n <- nrow(comp)
  ct <- tryCatch(stats::cor.test(comp$geneA, comp$geneB, method = method,
                                 exact = FALSE),
                 error = function(e) e)
  if (inherits(ct, "error")) {
    warning("viz_corplot(): cor.test failed: ", conditionMessage(ct),
            call. = FALSE)
    ct <- NULL
  }
  if (!is.null(ct) && !is.na(ct$estimate)) {
    rr <- round(as.numeric(ct$estimate), 3)
    pp <- round(ct$p.value, 3)
    pstr <- ifelse(pp == 0, "< 0.001", paste0("= ", pp))
    title_txt <- paste0("n = ", n, ", r = ", rr, " (", method,
                        "), p.value ", pstr)
  } else {
    title_txt <- paste0("n = ", n, ", r = NA (", method,
                        ") - correlation not estimable")
  }

  ggplot2::ggplot(comp, ggplot2::aes(geneA, geneB)) +
    ggplot2::geom_point(colour = "black") +
    ggplot2::geom_smooth(method = stats::lm, se = TRUE, na.rm = TRUE,
                         fullrange = TRUE, size = 2, colour = "red") +
    ggplot2::geom_rug(colour = "#006fbc") +
    ggplot2::theme_minimal() +
    ggplot2::xlab(paste0(a, x_lab)) +
    ggplot2::ylab(paste0(b, y_lab)) +
    ggplot2::labs(title = title_txt) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                   plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
                   axis.title = ggplot2::element_text(size = 16,
                                                      colour = "black"),
                   axis.text = ggplot2::element_text(size = 16,
                                                     colour = "black"))
}
