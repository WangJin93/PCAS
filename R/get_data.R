#' @title Get CPTAC data
#' @description
#' Get the CPTAC data by using the PCAS REST API
#' (\url{https://www.jingege.wang/bioinformatics/PCAS/api.php}).
#' This is the low-level entry point used by the other query functions.
#' @param table Table/dataset name. For \code{action = "expression"} use
#'   \code{dataset_info$Abbre}; for \code{action = "clinic"} remove the
#'   \code{_protein}/\code{_mRNA}/\code{_Phospho} suffix from
#'   \code{dataset_info$Abbre}.
#' @param action One of \code{"expression"}, \code{"DEGs"} or \code{"clinic"}.
#' @param genes Character vector of gene symbols / identifiers. Optional for
#'   \code{action = "DEGs"} and \code{action = "clinic"} (the server ignores it
#'   for those actions).
#' @param timeout Request timeout in seconds (default 60; large DEGs tables can
#'   occasionally be slow to assemble on the server).
#' @param tries Number of attempts before giving up (default 3).
#' @return A \code{data.frame} with the parsed API payload, or \code{NULL} when
#'   the request failed or no rows were returned (an explanatory message is
#'   always emitted first).
#' @details
#'   A failed or empty request returns \code{NULL} (never throws), so callers
#'   can uniformly test for missing data. Network/server errors are reported
#'   with \code{warning()}; an empty but successful response is reported with
#'   \code{message()}.
#' @examples
#' \dontrun{
#' results <- get_data(table = "LUAD_Academia_protein",
#'                     action = "expression",
#'                     genes = c("GAPDH", "TNS1"))
#' degs <- get_data(table = "LUAD_CPTAC_protein_limma", action = "DEGs")
#' clinic <- get_data(table = "LUAD_APOLLO", action = "clinic")
#' }
#' @export
get_data <- function(table = "LUAD_Academia_protein",
                     action = "expression",
                     genes = NULL,
                     timeout = 60,
                     tries = 3) {
  # ---------------------------------------------------------------------------
  # 1. Argument validation (early, with actionable feedback)
  # ---------------------------------------------------------------------------
  action <- match.arg(action, c("expression", "DEGs", "clinic"))

  if (!is.character(table) || length(table) != 1L || !nzchar(table)) {
    warning("get_data(): 'table' must be a single non-empty string.",
            call. = FALSE)
    return(NULL)
  }

  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    genes <- genes[nzchar(trimws(genes))]
    if (!length(genes)) genes <- NULL
  }
  if (action == "expression" && is.null(genes)) {
    warning("get_data(): 'genes' is required when action = \"expression\".",
            call. = FALSE)
    return(NULL)
  }

  timeout <- as.numeric(timeout)
  if (is.na(timeout) || timeout <= 0) timeout <- 60
  tries <- max(1L, as.integer(tries))

  # ---------------------------------------------------------------------------
  # 2. Build the URL (URL-encode each identifier, keep "," as the separator)
  # ---------------------------------------------------------------------------
  base <- "https://www.jingege.wang/bioinformatics/PCAS/api.php"
  q <- paste0("action=", utils::URLencode(action, reserved = TRUE),
              "&table=", utils::URLencode(table, reserved = TRUE))
  if (!is.null(genes)) {
    q <- paste0(q, "&genes=",
                paste(vapply(genes, utils::URLencode, character(1),
                             reserved = TRUE), collapse = ","))
  }
  url <- paste0(base, "?", q)

  # ---------------------------------------------------------------------------
  # 3. Request with timeout + retry
  # ---------------------------------------------------------------------------
  resp <- NULL
  last_err <- NULL
  for (i in seq_len(tries)) {
    ans <- tryCatch(
      httr::GET(url, httr::timeout(timeout),
                httr::user_agent("PCAS (R package)")),
      error = function(e) e
    )
    if (!inherits(ans, "error") && httr::status_code(ans) == 200L) {
      resp <- ans
      break
    }
    last_err <- ans
    if (i < tries) Sys.sleep(0.5 * i)
  }

  if (is.null(resp)) {
    detail <- if (inherits(last_err, "error")) {
      conditionMessage(last_err)
    } else if (!is.null(last_err)) {
      paste("HTTP", httr::status_code(last_err))
    } else {
      "unknown network error"
    }
    warning("get_data(): request failed for table '", table,
            "' (action = ", action, "): ", detail,
            ". Is the PCAS server reachable?", call. = FALSE)
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # 4. Parse JSON into a data.frame (server HTML / empty payloads handled)
  # ---------------------------------------------------------------------------
  txt <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"),
                  error = function(e) NULL)
  if (is.null(txt) || !nzchar(trimws(txt))) {
    warning("get_data(): empty response body from the PCAS server (table = '",
            table, "').", call. = FALSE)
    return(NULL)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(txt, simplifyDataFrame = TRUE),
    error = function(e) {
      warning("get_data(): could not parse the PCAS server response for table '",
              table, "' (", conditionMessage(e),
              "). The server may have returned an error page.", call. = FALSE)
      NULL
    }
  )
  if (is.null(parsed)) return(NULL)

  if (!is.data.frame(parsed)) {
    # e.g. an empty JSON array parsed to an empty list
    parsed <- NULL
  } else if (nrow(parsed) == 0L) {
    parsed <- NULL
  }

  if (is.null(parsed)) {
    .pcas_note("No rows returned for table '", table,
               "' (action = ", action, ").")
  }
  parsed
}
