#' @title Get CPTAC expression data
#' @description
#' Get the mRNA/protein/phosphosite expression data of one or several
#' identifiers in one or several CPTAC datasets.
#'
#' @param datasets Dataset names, you can input one or multiple datasets. Use
#'   \code{dataset_info$Abbre} to get all datasets.
#' @param genes Gene symbols / identifiers, one or multiple. For mRNA datasets,
#'   gene symbols are mapped internally through \code{idmap_RNA}; for
#'   phosphoproteome datasets pass phosphorylation-site ids (e.g.
#'   \code{"NP_000537.3:s315"}).
#' @param use_cache Logical; use the on-disk query cache (default TRUE).
#' @param cache_dir Optional custom cache directory. If NULL (default) the
#'   cache location is resolved from \code{options(PCAS.cache.dir)} and
#'   otherwise defaults to the per-user cache directory
#'   (\code{tools::R_user_dir("PCAS", "cache")}, i.e. \code{~/.cache/PCAS} on
#'   Linux). Set \code{use_cache = FALSE} or
#'   \code{options(PCAS.cache.dir = NA)} to disable caching entirely.
#' @return A data.frame with columns \code{ID}, \code{type}, \code{dataset} and
#'   one column per requested identifier (in the order of \code{genes}).
#'   Identifiers without measurements in a dataset are kept as all-NA columns so
#'   that downstream code never silently loses them. When no dataset returns any
#'   row, \code{NULL} is returned (with messages explaining why).
#' @details
#'   The return value carries an attribute \code{availability}, a
#'   \code{data.frame} with one row per \code{dataset x gene} combination giving
#'   \code{present} (measured or not) and \code{n_valid} (number of samples with
#'   a value). Read it with \code{attr(res, "availability")}.
#'
#'   If several mRNA probes map to the same gene symbol they are averaged per
#'   sample (a message is printed), so a sample never appears twice for the same
#'   gene.
#'
#'   \strong{Caching (same scheme as the GCAS package):} the processed result
#'   for every dataset is saved locally as
#'   \code{<cache_dir>/data_temp/<dataset>_<md5(ids)>.RData} and reused for the
#'   next identical request ("Loading cached data from ..."), so the PCAS server
#'   is not contacted again for the same data. Cache files do not expire
#'   automatically; delete the cache directory (or the individual file) to
#'   force a fresh download.
#' @examples
#' \dontrun{
#' results <- get_expr_data(datasets = "LUAD_CPTAC_mRNA",
#'                          genes = c("GAPDH", "TNS1"))
#' results <- get_expr_data(datasets = c("LUAD_CPTAC_protein",
#'                                        "LSCC_CPTAC_protein"),
#'                          genes = "GAPDH")
#' }
#' @export
get_expr_data <- function(datasets = c("LUAD_CPTAC_protein",
                                       "LSCC_CPTAC_protein"),
                          genes = c("TP53", "TNS1"),
                          use_cache = TRUE,
                          cache_dir = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Validate and normalise user input
  # ---------------------------------------------------------------------------
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) {
    .pcas_note("No gene identifiers supplied; returning NULL.")
    return(NULL)
  }

  datasets <- unique(trimws(as.character(datasets)))
  datasets <- datasets[nzchar(datasets)]
  if (!length(datasets)) {
    .pcas_note("No dataset names supplied; returning NULL.")
    return(NULL)
  }

  # dataset names must come from dataset_info (if it is visible)
  if (exists("dataset_info") && is.data.frame(dataset_info)) {
    unknown <- setdiff(datasets, dataset_info$Abbre)
    if (length(unknown)) {
      .pcas_note("Ignoring unknown dataset name(s) (use dataset_info$Abbre): ",
                 paste(unknown, collapse = ", "))
      datasets <- setdiff(datasets, unknown)
    }
    if (!length(datasets)) {
      .pcas_note("No valid dataset names remain; returning NULL.")
      return(NULL)
    }
  }

  cache_dir <- .pcas_resolve_cache_dir(use_cache, cache_dir)

  .pcas_note("Querying data of identifier ", paste(genes, collapse = ", "),
             " from datasets ", paste(datasets, collapse = ", "), ".")

  # ---------------------------------------------------------------------------
  # 2. Fetch dataset by dataset
  # ---------------------------------------------------------------------------
  parts   <- list()          # per-dataset long data.frame
  skipped <- character()     # datasets that produced nothing
  avail   <- list()          # availability records

  for (x in datasets) {
    is_mrna <- grepl("mRNA", x)
    map_df  <- NULL

    # -- mRNA datasets: map gene symbols to probe ids first ----------------
    if (is_mrna) {
      if (!exists("idmap_RNA") || !is.data.frame(idmap_RNA)) {
        stop("mRNA dataset '", x,
             "' requested but the idmap_RNA object is not available.",
             call. = FALSE)
      }
      map_df <- idmap_RNA[idmap_RNA$Symbol %in% genes,
                          c("row_names", "Symbol", "gene_type"), drop = FALSE]
      map_df <- map_df[!duplicated(map_df$row_names), , drop = FALSE]
      ids <- unique(map_df$row_names)
    } else {
      ids <- genes
    }

    record <- function(gene, present, n_valid, note) {
      avail[[length(avail) + 1L]] <<- data.frame(
        dataset = x, gene = gene, present = present, n_valid = n_valid,
        note = note, stringsAsFactors = FALSE)
    }

    # identifiers that can never be measured in this dataset
    if (is_mrna) {
      never <- setdiff(genes, map_df$Symbol)
      if (length(never)) {
        .pcas_note("Identifier(s) ", paste(never, collapse = ", "),
                   " were not found in the mRNA identifier map and are skipped ",
                   "for dataset ", x, ".")
        for (g in never) record(g, FALSE, 0L, "not in idmap_RNA")
      }
      if (!length(ids)) {
        skipped <- c(skipped, x)
        .pcas_note("No queryable identifier left for dataset ", x,
                   "; it is skipped.")
        next
      }
    }

    # -- cache (GCAS-style: <dataset>_<md5>.RData under <cache>/data_temp) ----
    cache_file <- NULL
    if (!is.null(cache_dir)) {
      cache_file <- .pcas_cache_file(cache_dir, x, ids, what = "data_temp")
      cache_file <- .pcas_ensure_cache_subdir(cache_file)
    }
    ds_df <- NULL
    if (!is.null(cache_file) && file.exists(cache_file)) {
      ds_df <- .pcas_load_cache(cache_file)
      if (!is.null(ds_df)) {
        .pcas_note("Loading cached data from ", cache_file, ".")
      }
    }

    # -- remote fetch --------------------------------------------------------
    if (is.null(ds_df)) {
      data <- get_data(x, "expression", ids)
      if (is.null(data)) {
        skipped <- c(skipped, x)
        .pcas_note("Dataset ", x, " returned no expression rows for the ",
                   "requested identifier(s); it is skipped.")
        for (g in genes) record(g, FALSE, 0L, "API returned no rows")
        next
      }

      # keep only the identifier column and the numeric sample columns
      keycol <- "row_names"
      num_df <- data[, setdiff(colnames(data), keycol), drop = FALSE]
      keys   <- as.character(data[[keycol]])

      # -- mRNA: attach symbols, collapse multi-probe symbols ----------------
      if (is_mrna) {
        mp <- map_df[, c("row_names", "Symbol"), drop = FALSE]
        keep <- keys %in% mp$row_names
        if (any(!keep)) {
          .pcas_note("Dropped ", sum(!keep),
                     " probe row(s) of dataset ", x,
                     " that had no entry in idmap_RNA.")
          keys   <- keys[keep]
          num_df <- num_df[keep, , drop = FALSE]
        }
        sym <- mp$Symbol[match(keys, mp$row_names)]
        keys <- sym
      }

      # collapse duplicated keys (multi-probe symbols / duplicated input)
      tab <- table(keys)
      multi <- names(tab)[tab > 1L]
      if (length(multi)) {
        .pcas_note("Averaging ", length(multi), " duplicated identifier row(s) ",
                   "(e.g. multiple mRNA probes for ",
                   paste(utils::head(multi, 5), collapse = ", "),
                   if (length(multi) > 5) ", ..." else "", ") in dataset ",
                   x, ".")
      }
      # transpose: build one value vector per sample per identifier
      sample_names <- colnames(num_df)
      smp <- .pcas_parse_sample(sample_names)
      ds_df <- smp
      for (g in genes) {
        idx <- which(keys == g)
        if (length(idx)) {
          mat <- suppressWarnings(
            vapply(idx, function(i) as.numeric(num_df[i, ]),
                   numeric(length(sample_names))))
          if (is.matrix(mat)) {
            v <- rowMeans(mat, na.rm = TRUE)
            v[is.nan(v)] <- NA_real_
          } else {
            v <- as.numeric(mat)
          }
          ds_df[[g]] <- v
        } else {
          ds_df[[g]] <- NA_real_
        }
      }
      rownames(ds_df) <- NULL
      ds_df$dataset <- x
      ds_df <- ds_df[, c("ID", "type", "dataset", genes), drop = FALSE]

      # write cache (best effort; a read-only FS must not abort the query)
      if (!is.null(cache_file)) .pcas_save_cache(ds_df, cache_file)
    }

    parts[[x]] <- ds_df
    for (g in genes) {
      record(g, any(!is.na(ds_df[[g]])), sum(!is.na(ds_df[[g]])),
             if (is_mrna && !(g %in% map_df$Symbol)) "not in idmap_RNA" else "")
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Combine and report
  # ---------------------------------------------------------------------------
  if (!length(parts)) {
    .pcas_note("No expression data could be retrieved for any dataset; ",
               "returning NULL.")
    return(NULL)
  }
  combined <- dplyr::bind_rows(parts)

  # keep only sample rows carrying at least one measured value
  n_before <- nrow(combined)
  gene_cols <- intersect(genes, colnames(combined))
  has_val <- rowSums(!is.na(combined[, gene_cols, drop = FALSE])) > 0L
  combined <- combined[has_val, , drop = FALSE]
  if (nrow(combined) < n_before) {
    .pcas_note("Dropped ", n_before - nrow(combined),
               " sample row(s) with no value for any requested identifier.")
  }
  if (!nrow(combined)) {
    .pcas_note("No sample carried any value for the requested identifier(s); ",
               "returning NULL.")
    return(NULL)
  }

  combined <- combined[, c("ID", "type", "dataset", genes), drop = FALSE]
  for (g in genes) {
    combined[[g]] <- .pcas_as_numeric(combined[[g]])
  }
  rownames(combined) <- NULL

  availability <- dplyr::bind_rows(avail)
  missing_all <- genes[vapply(genes, function(g) {
    !any(!is.na(combined[[g]]))
  }, logical(1))]
  if (length(missing_all)) {
    .pcas_note("Identifier(s) with no measurement in any dataset: ",
               paste(missing_all, collapse = ", "))
  }
  if (length(skipped)) {
    .pcas_note("Dataset(s) skipped because they returned no data: ",
               paste(skipped, collapse = ", "))
  }
  attr(combined, "availability") <- availability
  .pcas_note("Returned expression for ", nrow(combined), " sample(s) across ",
             length(unique(combined$dataset)), " dataset(s).")
  combined
}
