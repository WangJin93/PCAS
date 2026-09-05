#' @title Get DEGs results
#' @description
#' Get the results of differential expression analysis between tumour and
#' normal samples in a CPTAC dataset (pre-computed with limma or t-test on the
#' PCAS server).
#' @param dataset Dataset abbreviation, e.g. \code{"LUAD_CPTAC_protein"}. Use
#'   \code{dataset_info$Abbre} to list all datasets.
#' @param method One of \code{"t.test"} or \code{"limma"}.
#' @param use_cache Logical; use the on-disk query cache (default TRUE).
#' @return A data.frame of differential expression results with at least the
#'   columns \code{Symbol}, \code{logFC}, \code{P.Value} and
#'   \code{adj.P.Val}, or \code{NULL} when the dataset is unknown or the server
#'   returned no rows (messages/warnings explain why). The returned object
#'   carries \code{attr(x, "dataset")} and \code{attr(x, "method")}.
#' @examples
#' \dontrun{
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "limma")
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_mRNA", method = "t.test")
#' }
#' @export
get_DEGs_result <- function(dataset = "LUAD_CPTAC_protein",
                            method = "t.test",
                            use_cache = TRUE) {
  method <- match.arg(method, c("t.test", "limma"))

  if (!is.character(dataset) || length(dataset) != 1L || !nzchar(dataset)) {
    warning("get_DEGs_result(): 'dataset' must be a single non-empty string.",
            call. = FALSE)
    return(NULL)
  }

  # dataset must be a real expression table from dataset_info
  if (exists("dataset_info") && is.data.frame(dataset_info)) {
    if (!dataset %in% dataset_info$Abbre) {
      warning("get_DEGs_result(): dataset '", dataset,
              "' is not a known CPTAC table. Use dataset_info$Abbre to list ",
              "valid datasets (e.g. 'LUAD_CPTAC_protein').", call. = FALSE)
      return(NULL)
    }
  }

  table <- paste0(dataset, ifelse(method == "limma", "_limma", "_ttest"))

  cache_file <- NULL
  cache_dir  <- .pcas_cache_dir(use_cache)
  if (!is.null(cache_dir)) {
    cache_file <- .pcas_cache_file(cache_dir, table, "DEGs", what = "DEGs")
    cache_file <- .pcas_ensure_cache_subdir(cache_file)
  }

  results <- NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {
    e <- new.env(parent = emptyenv())
    ok <- tryCatch({load(cache_file, envir = e); TRUE},
                   error = function(e) FALSE)
    if (ok && exists("results", envir = e, inherits = FALSE)) {
      results <- e$results
      .pcas_note("Loaded cached DEGs results for ", table, ".")
    }
  }

  if (is.null(results)) {
    .pcas_note("Querying DEGs results for table ", table, ".")
    results <- get_data(table, action = "DEGs")
    if (is.null(results)) {
      warning("get_DEGs_result(): no DEGs data returned for dataset '",
              dataset, "' (method ", method, ").", call. = FALSE)
      return(NULL)
    }
    n_before <- nrow(results)

    # -- mRNA datasets: first column holds probe ids; attach gene symbols -----
    if (grepl("mRNA", dataset)) {
      colnames(results)[1] <- "mRNAs"
      # The server returns probe ids that may be longer than the 15-character
      # ids stored in idmap_RNA; truncate to stay mergeable.
      results$mRNAs <- substr(results$mRNAs, 1, 15)
      id_map <- idmap_RNA[, c("mRNAs", "Symbol", "gene_type"), drop = FALSE]
      id_map <- id_map[!duplicated(id_map$mRNAs), , drop = FALSE]
      merged <- merge(id_map, results, by = "mRNAs")
      if (nrow(merged) < n_before) {
        .pcas_note("mRNA DEGs: dropped ", n_before - nrow(merged),
                   " row(s) whose probe id could not be mapped to a gene ",
                   "symbol in idmap_RNA.")
      }
      results <- merged
      if (!nrow(results)) {
        warning("get_DEGs_result(): no DEGs row could be mapped to a gene ",
                "symbol for dataset '", dataset, "'.", call. = FALSE)
        return(NULL)
      }
    } else {
      colnames(results)[1] <- "Symbol"
    }

    # -- numeric coercion with feedback ---------------------------------------
    for (col in intersect(c("logFC", "P.Value", "adj.P.Val"),
                          colnames(results))) {
      results[[col]] <- .pcas_as_numeric(results[[col]])
    }
    if (all(is.na(results$logFC)) || all(is.na(results$P.Value))) {
      warning("get_DEGs_result(): the returned DEGs table for '", dataset,
              "' has no usable logFC/P.Value values.", call. = FALSE)
      return(NULL)
    }

    if (!is.null(cache_file)) {
      tryCatch(save(results, file = cache_file), error = function(e) {
        .pcas_note("Could not write cache file ", cache_file,
                   " (", conditionMessage(e), "); continuing without it.")
      })
    }
  }

  attr(results, "dataset") <- dataset
  attr(results, "method")  <- method
  results
}
