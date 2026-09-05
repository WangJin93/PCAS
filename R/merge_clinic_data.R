#' @title Merge clinic data
#' @description
#' Fetch the clinical table of a cohort and merge it with expression data
#' (tumour samples only). Contrary to a silent inner join, the function reports
#' exactly how many samples matched and how many clinical fields are missing.
#'
#' @param cohort Data cohort, for example \code{"LUAD_APOLLO"}, \code{"LUAD_CPTAC"}
#'   (a trailing \code{_protein}/\code{_mRNA}/\code{_Phospho} suffix, if
#'   present, is stripped automatically).
#' @param data_input Expression data obtained from \code{\link{get_expr_data}()}.
#' @param return_summary Logical; if TRUE (default) a list
#'   \code{list(df = <merged data.frame>, summary = <matching/NA summary>)} is
#'   returned. If FALSE the plain merged data.frame is returned with the summary
#'   attached as \code{attr(x, "merge_summary")} for backward compatibility.
#' @return Either a list with the merged data.frame and a summary, or (when
#'   \code{return_summary = FALSE}) the merged data.frame alone. Returns
#'   \code{NULL} with explanatory messages when the clinical table cannot be
#'   fetched or nothing can be merged.
#' @details
#'   The summary contains: input sample counts, matched/unmatched sample counts,
#'   per-clinical-field missingness, and which simplified stage columns were
#'   added. Rows are never silently dropped without being accounted for in the
#'   summary.
#' @examples
#' \dontrun{
#' data_input <- get_expr_data("LUAD_APOLLO_mRNA", "TP53")
#' res <- merge_clinic_data("LUAD_APOLLO", data_input)
#' head(res$df)          # merged expression + clinical columns
#' res$summary           # matching / missingness report
#' }
#' @export
merge_clinic_data <- function(cohort = "LUAD_APOLLO",
                              data_input,
                              return_summary = TRUE) {
  # ---------------------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data_input)) {
    stop("merge_clinic_data(): 'data_input' must be a data.frame obtained from ",
         "get_expr_data(). Got: ", paste(class(data_input), collapse = "/"),
         call. = FALSE)
  }
  .pcas_require_cols(data_input, c("ID", "type"), "merge_clinic_data")

  cohort <- as.character(cohort)
  cohort <- .pcas_strip_dataset(cohort[1L])

  # ---------------------------------------------------------------------------
  # 2. Fetch clinical data
  # ---------------------------------------------------------------------------
  clinic <- get_data(cohort, "clinic")
  if (is.null(clinic)) {
    warning("merge_clinic_data(): no clinical data returned for cohort '",
            cohort, "'. Check the cohort name (e.g. 'LUAD_APOLLO').",
            call. = FALSE)
    return(NULL)
  }
  if (!"Cases_Submitter_ID" %in% colnames(clinic)) {
    stop("merge_clinic_data(): the clinical table of cohort '", cohort,
         "' has no 'Cases_Submitter_ID' column.", call. = FALSE)
  }

  # drop the generic row_name column if present
  clinic <- clinic[, setdiff(colnames(clinic), "row_names"), drop = FALSE]

  # ---------------------------------------------------------------------------
  # 3. Filter tumour samples and merge with clinical records
  # ---------------------------------------------------------------------------
  n_expr_rows  <- nrow(data_input)
  n_expr_id    <- length(unique(data_input$ID))
  tumour       <- data_input[data_input$type == "Tumor", , drop = FALSE]
  n_tumour_rows <- nrow(tumour)
  n_tumour_id  <- length(unique(tumour$ID))
  if (!n_tumour_rows) {
    .pcas_note("data_input contains no 'Tumor' sample; nothing to merge for ",
               "cohort ", cohort, ".")
  }

  n_clinic <- length(unique(clinic$Cases_Submitter_ID))
  merged   <- merge(tumour, clinic, by.x = "ID",
                    by.y = "Cases_Submitter_ID", all.x = FALSE)

  n_match_rows <- nrow(merged)
  n_unmatch    <- n_tumour_rows - n_match_rows

  # ---------------------------------------------------------------------------
  # 4. Simplified stage columns (only when the source column exists)
  # ---------------------------------------------------------------------------
  simplify_map <- c(
    "AJCC_Pathologic_Stage" = "AJCC_Pathologic_Stage_simplify",
    "Tumor_Stage"           = "Tumor_Stage_simplify",
    "AJCC_Pathologic_T"     = "AJCC_Pathologic_T_simplify"
  )
  added <- character(0)
  for (src in names(simplify_map)) {
    if (src %in% colnames(merged)) {
      dst <- simplify_map[[src]]
      pat <- if (grepl("AJCC_Pathologic_T", src, fixed = TRUE)) "[abc]" else "[ABC]"
      merged[[dst]] <- sub(pat, "", merged[[src]])
      added <- c(added, dst)
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Missingness report for the merged clinical fields
  # ---------------------------------------------------------------------------
  clin_cols <- setdiff(colnames(clinic), "Cases_Submitter_ID")
  na_report <- data.frame(
    field     = clin_cols,
    n_missing = vapply(clin_cols, function(cl) sum(is.na(merged[[cl]])),
                       integer(1)),
    pct_missing = vapply(clin_cols, function(cl) {
      100 * mean(is.na(merged[[cl]]))
    }, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  summary <- list(
    cohort                = cohort,
    n_expr_samples        = n_expr_rows,
    n_expr_unique_ids     = n_expr_id,
    n_tumour_rows         = n_tumour_rows,
    n_tumour_unique_ids   = n_tumour_id,
    n_clinic_patients     = n_clinic,
    n_matched_rows        = n_match_rows,
    n_unmatched_tumour_rows = n_unmatch,
    clinical_na           = na_report,
    simplified_stage_cols = added
  )

  # ---------------------------------------------------------------------------
  # 6. Feedback + return
  # ---------------------------------------------------------------------------
  if (n_unmatch > 0L) {
    .pcas_note(n_unmatch, " tumour expression row(s) (of ", n_tumour_rows,
               ") of cohort ", cohort,
               " have no clinical record and were left out of the merge.")
  }
  high_na <- na_report$field[na_report$pct_missing >= 50 &
                               na_report$n_missing > 0]
  if (length(high_na)) {
    .pcas_note("Clinical field(s) missing in >= 50% of matched samples: ",
               paste(high_na, collapse = ", "), ".")
  }
  if (length(added)) {
    .pcas_note("Added simplified stage column(s): ",
               paste(added, collapse = ", "), ".")
  }

  if (return_summary) {
    list(df = merged, summary = summary)
  } else {
    attr(merged, "merge_summary") <- summary
    merged
  }
}
