#' @title Get DEGs results
#' @description
#' Get the results of different expression analysis between tumor and normal samples in CPTAC datasets.
#' @param dataset Use dataset$Abbre to get all tables.
#' @param method "limma" or "t.test".
#' @examples
#' \dontrun{
#' results <- get_DEGs_result(dataset = "LUAD_CPTAC_protein", method = "limma")
#' }
#' @export
#'
get_DEGs_result <- function(dataset = "LUAD_CPTAC_protein",
                            method = "t.test"){
  table <- ifelse(method == "limma",paste0(dataset,"_limma"),paste0(dataset,"_ttest"))
  results <- get_data( table, action = "DEGs")
  if (stringr::str_detect(dataset,"mRNA")){
    colnames(results)[1] <- "mRNAs"
    results$mRNAs <- substr(results$mRNAs,1,15)
    results <- merge(idmap_RNA[c(1,3,4)],results,by="mRNAs")
  }else{
    colnames(results)[1] <- "Symbol"
  }
  results$logFC <- as.numeric(results$logFC)
  results$P.Value <- as.numeric(results$P.Value)

  return(results)
}
