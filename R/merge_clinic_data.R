#' @title Merge clinic data
#' @description
#' Get clinic data and merge it with expression data.
#' @import dplyr
#' @param cohort Data cohort, for example, "LUAD_APOLLO", "LUAD_CPTAC".
#' @param data_input Expression data obtained from get_expr_data() function.
#' @examples
#' \dontrun{
#' data_input <- get_expr_data("LUAD_APOLLO_mRNA","TP53")
#' results <- merge_clinic_data("LUAD_APOLLO",data_input)
#' }
#' @export
#'
merge_clinic_data <-function(cohort="LUAD_APOLLO",data_input){
  clinic <- get_data(cohort,
                     "clinic")
  merge_data <- data_input  %>%
    dplyr::filter(type == "Tumor")  %>%
    merge(.,clinic[-1],by.x = "ID", by.y = "Cases_Submitter_ID")
  if ("AJCC_Pathologic_Stage" %in% colnames(merge_data)) merge_data$AJCC_Pathologic_Stage_simplify <-  stringr::str_remove( merge_data[,"AJCC_Pathologic_Stage"] ,"A|B|C")
  if ("Tumor_Stage" %in% colnames(merge_data)) merge_data$Tumor_Stage_simplify <- stringr::str_remove( merge_data[,"Tumor_Stage"] ,"A|B|C")
  if ("AJCC_Pathologic_T" %in% colnames(merge_data)) merge_data$AJCC_Pathologic_T_simplify <-  stringr::str_remove( merge_data[,"AJCC_Pathologic_T"] ,"a|b|c")
  return(merge_data)
}
