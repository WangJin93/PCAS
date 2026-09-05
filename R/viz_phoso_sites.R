#' @title Phosphorylation site visualisation on a protein diagram
#' @description
#' Draw the domain structure of a protein (via the UniProt features API,
#' through \pkg{drawProteins}) and overlay its phosphorylation sites. Sites are
#' taken either from the CPTAC phosphoproteomics table (\code{idmap_protein},
#' default) or from the UniProt features themselves.
#' @param gene Gene/protein symbol, e.g. \code{"TNS1"}.
#' @param phoso_infoDB One of \code{"CPTAC"} (default) or \code{"UniProt"}.
#' @return A ggplot object. The plotted site table (if any) is attached as
#'   \code{attr(p, "sites")} for export/verification. \code{NULL} (with an
#'   explanatory message) is returned when the gene cannot be mapped or the
#'   UniProt network call fails.
#' @details
#'   Combined CPTAC site identifiers such as \code{NP_000025.1:s218y223t227}
#'   are parsed into all of their sites (s218, y223, t227), not only the last
#'   one.
#' @examples
#' \dontrun{
#' viz_phoso_sites("TNS1")
#' viz_phoso_sites("YTHDC2", phoso_infoDB = "UniProt")
#' }
#' @export
viz_phoso_sites <- function(gene = "YTHDC2",
                            phoso_infoDB = "CPTAC") {
  phoso_infoDB <- match.arg(phoso_infoDB, c("UniProt", "CPTAC"))
  gene <- trimws(as.character(gene)[1L])
  if (is.na(gene) || !nzchar(gene)) {
    .pcas_note("viz_phoso_sites(): 'gene' must be a non-empty symbol.")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # 1. Map gene symbol -> reviewed UniProt entry
  # ---------------------------------------------------------------------------
  if (!exists("uniport_map") || !is.data.frame(uniport_map)) {
    stop("viz_phoso_sites(): the uniport_map object is not available.",
         call. = FALSE)
  }
  hits <- uniport_map[uniport_map$Symbol == gene, , drop = FALSE]
  reviewed <- unique(hits$Entry[hits$Reviewed == "reviewed"])
  if (!length(reviewed)) {
    .pcas_note("viz_phoso_sites(): no reviewed UniProt entry found for gene '",
               gene, "'; nothing can be drawn.")
    return(NULL)
  }
  uni_id <- reviewed[1L]
  if (length(reviewed) > 1L) {
    .pcas_note("Gene '", gene, "' maps to several reviewed UniProt entries; ",
               "drawing the first one (", uni_id, ").")
  }

  # ---------------------------------------------------------------------------
  # 2. Fetch UniProt features (network) and build the protein canvas
  # ---------------------------------------------------------------------------
  rel_json <- tryCatch(drawProteins::get_features(uni_id),
                       error = function(e) e)
  if (inherits(rel_json, "error") || is.null(rel_json) ||
      !length(rel_json) || length(rel_json) == 1L &&
      all(nchar(trimws(unlist(rel_json))) == 0L)) {
    warning("viz_phoso_sites(): could not fetch UniProt features for ", uni_id,
            " (", if (inherits(rel_json, "error")) conditionMessage(rel_json)
              else "empty response",
            "). Check your internet connection to UniProt.", call. = FALSE)
    return(NULL)
  }
  rel_data <- tryCatch(drawProteins::feature_to_dataframe(rel_json),
                       error = function(e) e)
  if (inherits(rel_data, "error")) {
    warning("viz_phoso_sites(): feature_to_dataframe failed: ",
            conditionMessage(rel_data), call. = FALSE)
    return(NULL)
  }
  if (!nrow(rel_data)) {
    .pcas_note("The UniProt record of '", gene,
               "' contains no structured feature; nothing can be drawn.")
    return(NULL)
  }

  p <- drawProteins::draw_canvas(rel_data)
  p <- drawProteins::draw_chains(p, rel_data, label_size = 5, labels = gene)

  # ---------------------------------------------------------------------------
  # 3. Phosphorylation sites
  # ---------------------------------------------------------------------------
  site_df <- NULL
  if (phoso_infoDB == "CPTAC") {
    if (!exists("idmap_protein") || !is.data.frame(idmap_protein)) {
      stop("viz_phoso_sites(): the idmap_protein object is not available.",
           call. = FALSE)
    }
    phoso_data <- idmap_protein[idmap_protein$Symbol == gene, , drop = FALSE]
    ids <- phoso_data$row_names
    ids <- ids[grepl(":", ids, fixed = TRUE)]        # keep phosphosite rows only

    records <- list()
    for (id in ids) {
      prefix    <- sub(":.*$", "", id)
      site_part <- sub("^[^:]*:", "", id)
      tokens    <- regmatches(site_part,
                              gregexpr("[sSyYtT][0-9]+", site_part))[[1L]]
      if (!length(tokens)) next
      records[[length(records) + 1L]] <- data.frame(
        aa        = toupper(substr(tokens, 1L, 1L)),
        phoso_site = paste0(prefix, ":", tokens),
        location  = suppressWarnings(as.numeric(substring(tokens, 2L))),
        order     = 1L,
        stringsAsFactors = FALSE)
    }
    if (length(records)) {
      site_df <- unique(do.call(rbind, records))
      .pcas_note("Found ", nrow(site_df), " CPTAC phosphorylation site(s) ",
                 "for gene '", gene, "'.")
    } else {
      .pcas_note("No CPTAC phosphorylation site found for gene '", gene,
                 "'; drawing the protein structure only.")
    }
  } else {
    # UniProt mode: look for phospho/modification features in rel_data
    if (is.data.frame(rel_data) && nrow(rel_data)) {
      type_col <- if ("type" %in% colnames(rel_data)) "type" else NA_character_
      if (!is.na(type_col)) {
        is_phos <- grepl("PHOS", rel_data[[type_col]], ignore.case = TRUE)
        if ("description" %in% colnames(rel_data)) {
          is_phos <- is_phos |
            (rel_data[[type_col]] == "MOD_RES" &
               grepl("[Pp]hospho", rel_data$description))
        }
        cand <- rel_data[is_phos, , drop = FALSE]
        if (nrow(cand) && all(c("begin", "order") %in% colnames(cand))) {
          site_df <- data.frame(aa = NA_character_,
                                phoso_site = cand$type,
                                location = cand$begin,
                                order = cand$order,
                                stringsAsFactors = FALSE)
        }
      }
    }
    if (is.null(site_df)) {
      .pcas_note("No phosphorylated residue feature found in the UniProt ",
                 "record of '", gene, "'; drawing the protein structure only.")
    }
  }

  if (!is.null(site_df) && nrow(site_df)) {
    site_df <- site_df[is.finite(site_df$location), , drop = FALSE]
    p <- p +
      ggplot2::geom_segment(data = site_df,
                            ggplot2::aes(x = location, y = order + 0.2,
                                         xend = location, yend = order + 0.3),
                            colour = "grey50", linewidth = 0.5, linetype = 1) +
      ggplot2::geom_point(data = site_df,
                          ggplot2::aes(x = location, y = order + 0.3),
                          shape = 21, colour = "black", fill = "red",
                          size = 4, show.legend = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 4. Draw remaining structural layers + cosmetics
  # ---------------------------------------------------------------------------
  p <- tryCatch(drawProteins::draw_domains(p, rel_data, label_domains = FALSE),
                error = function(e) {
                  .pcas_note("draw_domains skipped: ", conditionMessage(e))
                  p
                })
  p <- tryCatch(drawProteins::draw_regions(p, rel_data), error = function(e) {
    .pcas_note("draw_regions skipped: ", conditionMessage(e))
    p
  })
  p <- tryCatch(drawProteins::draw_motif(p, rel_data), error = function(e) {
    .pcas_note("draw_motif skipped: ", conditionMessage(e))
    p
  })

  p <- p + ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "bottom")

  if (!is.null(site_df)) attr(p, "sites") <- site_df
  p
}
