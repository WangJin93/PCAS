#' @title Visualize gene expression (tumour vs normal)
#' @description
#' Visualize expression differences between Tumour and Normal samples of CPTAC
#' data with box/violin plots: one gene of one dataset (\code{"single"}), several
#' genes of one dataset (\code{"multi_gene"}), or one gene across several
#' datasets (\code{"multi_set"}).
#' @param df Gene expression data obtained from \code{\link{get_expr_data}()}.
#' @param df_type One of \code{"single"}, \code{"multi_gene"},
#'   \code{"multi_set"}.
#' @param Show.P.value Whether to display the result of the differential test,
#'   default TRUE.
#' @param Show.P.label Whether to display significance markers
#'   ("***" / "**" / "*" / "ns") instead of numeric p-values, default TRUE.
#' @param Method Method of the differential analysis, \code{"t.test"} or
#'   \code{"wilcox.test"} (anything accepted by
#'   \code{ggpubr::compare_means}).
#' @param values Color palette for the two groups, default
#'   \code{c("#00AFBB", "#FC4E07")} (tumour, normal).
#' @param Show.n Display sample sizes below the plot.
#' @param Show.n.location Y position used for the sample-size labels; default
#'   \code{"default"} places them below the data range.
#' @return A ggplot object. Whenever rows are removed or a comparison is
#'   impossible (e.g. a dataset has no Normal samples), messages explain what
#'   happened instead of failing silently.
#' @examples
#' \dontrun{
#' df_single <- get_expr_data(datasets = "LUAD_CPTAC_mRNA", genes = "TP53")
#' df_multi_gene <- get_expr_data(datasets = "LUAD_CPTAC_protein",
#'                                genes = c("TP53", "TNS1"))
#' viz_TvsN(df_single, df_type = "single")
#' viz_TvsN(df_multi_gene, df_type = "multi_gene")
#' }
#' @export
viz_TvsN <- function(df,
                     df_type = "single",
                     Show.P.value = TRUE,
                     Show.P.label = TRUE,
                     Method = "t.test",
                     values = c("#00AFBB", "#FC4E07"),
                     Show.n = TRUE,
                     Show.n.location = "default") {
  # ---------------------------------------------------------------------------
  # validation
  # ---------------------------------------------------------------------------
  df_type <- match.arg(df_type, c("single", "multi_gene", "multi_set"))
  if (is.null(df) || !is.data.frame(df)) {
    .pcas_note("viz_TvsN(): 'df' is NULL/not a data.frame; returning NULL.")
    return(NULL)
  }
  .pcas_require_cols(df, c("ID", "type", "dataset"), "viz_TvsN")
  meta <- c("ID", "type", "dataset")
  expr_cols <- setdiff(colnames(df), meta)
  if (!length(expr_cols)) {
    .pcas_note("viz_TvsN(): 'df' contains no expression column.")
    return(NULL)
  }
  Method <- Method[1L]

  # keep only the two comparison groups, with feedback
  n_dropped_other <- sum(!df$type %in% c("Tumor", "Normal"))
  if (n_dropped_other > 0L) {
    .pcas_note("Removed ", n_dropped_other, " sample row(s) whose type is ",
               "neither Tumour nor Normal (types found: ",
               paste(unique(df$type), collapse = ", "), ").")
  }
  df <- df[df$type %in% c("Tumor", "Normal"), , drop = FALSE]
  if (!nrow(df)) {
    .pcas_note("viz_TvsN(): no Tumour/Normal sample row left.")
    return(NULL)
  }

  run_pvalue <- function(sub, grp) {
    ## safe compare_means wrapper: returns data.frame or NULL
    lv <- unique(sub$type)
    tb <- table(sub$type)
    if (length(lv) != 2L || any(tb < 2L) || all(is.na(sub[["value"]]))) {
      .pcas_note("Cannot test '", grp,
                 "': need >= 2 samples of each of Tumour and Normal; ",
                 "p-value is skipped for it.")
      return(NULL)
    }
    tryCatch(
      ggpubr::compare_means(value ~ type, data = sub, method = Method),
      error = function(e) {
        .pcas_note("Differential test failed for '", grp, "': ",
                   conditionMessage(e), "; p-value skipped.")
        NULL
      }
    )
  }

  ymax_global <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) 1 else max(v)
  }
  ymin_global <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) 0 else min(v)
  }

  # ===========================================================================
  if (df_type == "single") {
    value_col <- expr_cols[1L]
    if (length(expr_cols) > 1L) {
      .pcas_note("df contains ", length(expr_cols),
                 " expression columns; plotting only '", value_col, "'.")
    }
    df$value <- .pcas_as_numeric(df[[value_col]])
    df$gene  <- value_col

    na_rows <- is.na(df$value)
    if (any(na_rows)) {
      .pcas_note("Removed ", sum(na_rows), " row(s) with missing expression ",
                 "of '", value_col, "'.")
      df <- df[!na_rows, , drop = FALSE]
    }
    if (!nrow(df)) {
      .pcas_note("viz_TvsN(): no valid expression value left.")
      return(NULL)
    }
    pv <- if (Show.P.value) run_pvalue(df, value_col) else NULL

    count_N <- as.data.frame(table(df$type), stringsAsFactors = FALSE)
    names(count_N) <- c("type", "n")
    count_N$label <- paste0("n = ", count_N$n)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = type, y = value, fill = type)) +
      ggplot2::geom_violin(trim = FALSE, show.legend = FALSE) +
      ggplot2::geom_boxplot(width = 0.2, fill = "white", show.legend = FALSE) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(paste0(value_col, " expression")) +
      ggplot2::scale_fill_manual(values = values)

    ymin <- ymin_global(df$value); ymax <- ymax_global(df$value)
    if (Show.n) {
      nloc <- if (identical(Show.n.location, "default")) {
        ymin - (ymax - ymin) * 0.2
      } else Show.n.location
      p <- p +
        ggplot2::geom_text(data = count_N,
                           ggplot2::aes(x = type, y = nloc, label = label),
                           colour = values[seq_len(nrow(count_N))],
                           size = 6, hjust = 0.5, show.legend = FALSE)
    }
    if (Show.P.value && !is.null(pv) && nrow(pv)) {
      yp <- ymax + (ymax - ymin) * 0.1
      lab <- if (Show.P.label) {
        pv$p.signif
      } else {
        ifelse(signif(pv$p, 3) < 0.001, "P < 0.001",
               paste0("P = ", signif(pv$p, 3)))
      }
      p <- p + ggplot2::annotate("text", x = 1.5, y = yp, label = lab,
                                 size = if (Show.P.label) 8 else 6)
    }
  }

  # ===========================================================================
  if (df_type == "multi_gene") {
    long <- reshape2::melt(df[, c(meta, expr_cols), drop = FALSE],
                           id.vars = meta,
                           measure.vars = expr_cols,
                           variable.name = "gene",
                           value.name = "value")

    long$value <- .pcas_as_numeric(long$value)
    na_rows <- is.na(long$value)
    if (any(na_rows)) {
      .pcas_note("Removed ", sum(na_rows), " gene-sample row(s) with missing ",
                 "expression before plotting.")
      long <- long[!na_rows, , drop = FALSE]
    }
    if (!nrow(long)) {
      .pcas_note("viz_TvsN(): no valid expression value left.")
      return(NULL)
    }

    pv <- NULL
    if (Show.P.value) {
      parts <- lapply(unique(long$gene), function(g) {
        res <- run_pvalue(long[long$gene == g, , drop = FALSE], g)
        if (!is.null(res)) {
          res$gene <- g
          res
        }
      })
      pv <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
    }

    count_N <- as.data.frame(table(long$gene, long$type),
                             stringsAsFactors = FALSE)
    names(count_N) <- c("gene", "type", "n")
    count_N$label <- paste0("n = ", count_N$n)

    p <- ggplot2::ggplot(long,
                         ggplot2::aes(x = gene, y = value, fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("Expression") +
      ggplot2::scale_fill_manual(values = values)

    ymin <- ymin_global(long$value); ymax <- ymax_global(long$value)
    if (Show.n) {
      nloc <- if (identical(Show.n.location, "default")) {
        ymin - (ymax - ymin) * 0.2
      } else Show.n.location
      p <- p + ggplot2::geom_text(
        data = count_N,
        ggplot2::aes(x = gene, y = nloc, label = label, colour = type),
        position = ggplot2::position_dodge2(0.9),
        size = 6, angle = 90, hjust = 0, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = values)
    }
    if (Show.P.value && !is.null(pv) && nrow(pv)) {
      yp <- ymax + (ymax - ymin) * 0.1
      lab <- if (Show.P.label) {
        pv$p.signif
      } else {
        paste0("P = ", signif(pv$p, 2))
      }
      p <- p + ggplot2::geom_text(
        data = pv,
        ggplot2::aes(x = gene, y = yp, label = lab), inherit.aes = FALSE,
        size = if (Show.P.label) 6 else 5)
    }
  }

  # ===========================================================================
  if (df_type == "multi_set") {
    value_col <- expr_cols[1L]
    if (length(expr_cols) > 1L) {
      .pcas_note("df contains ", length(expr_cols),
                 " expression columns; plotting only '", value_col, "'.")
    }
    df$value <- .pcas_as_numeric(df[[value_col]])
    na_rows <- is.na(df$value)
    if (any(na_rows)) {
      .pcas_note("Removed ", sum(na_rows), " row(s) with missing expression ",
                 "of '", value_col, "'.")
      df <- df[!na_rows, , drop = FALSE]
    }
    if (!nrow(df)) {
      .pcas_note("viz_TvsN(): no valid expression value left.")
      return(NULL)
    }

    pv <- NULL
    if (Show.P.value) {
      parts <- lapply(unique(df$dataset), function(ds) {
        res <- run_pvalue(df[df$dataset == ds, , drop = FALSE], ds)
        if (!is.null(res)) {
          res$dataset <- ds
          res
        }
      })
      pv <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
    }

    count_N <- as.data.frame(table(df$dataset, df$type),
                             stringsAsFactors = FALSE)
    names(count_N) <- c("dataset", "type", "n")
    count_N$label <- paste0("n = ", count_N$n)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = dataset, y = value,
                                          fill = type)) +
      ggplot2::geom_boxplot() +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(paste0(value_col, " expression")) +
      ggplot2::scale_fill_manual(values = values)

    ymin <- ymin_global(df$value); ymax <- ymax_global(df$value)
    if (Show.n) {
      nloc <- if (identical(Show.n.location, "default")) {
        ymin - (ymax - ymin) * 0.2
      } else Show.n.location
      p <- p + ggplot2::geom_text(
        data = count_N,
        ggplot2::aes(x = dataset, y = nloc, label = label, colour = type),
        position = ggplot2::position_dodge2(0.9),
        size = 5, angle = 90, hjust = 0, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = values)
    }
    if (Show.P.value && !is.null(pv) && nrow(pv)) {
      yp <- ymax * 1.1
      lab <- if (Show.P.label) {
        pv$p.signif
      } else {
        paste0("P = ", signif(pv$p, 2))
      }
      p <- p + ggplot2::geom_text(
        data = pv, ggplot2::aes(x = dataset, y = yp, label = lab),
        inherit.aes = FALSE, size = if (Show.P.label) 6 else 5)
    }
  }

  p + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 16),
                   axis.title.y = ggplot2::element_text(size = 16)) +
    (if (df_type == "multi_set") {
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1.0),
        plot.margin = ggplot2::unit(c(0.2, 0.2, 0.2, 2), "cm"))
    } else ggplot2::theme())
}
