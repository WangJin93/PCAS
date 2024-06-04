#' @title Get CPTAC data
#' @description
#' Get the CPTAC data by using the api. All results saved in MySQL database.
#' @import jsonlite
#' @param table For action = expression, use dataset$Abbre to get all tables; For action = clinic, remove _protein/_mRNA/_Phospho from dataset$Abbre.
#' @param action "expression", "degs" or "clinic".
#' @param gene Gene symbols, you can input one or multiple symbols.
#' @examples
#' \dontrun{
#' results <- get_data(table = "LUAD_Academia_protein", action = "expression", genes = c("GAPDH","TNS1"))
#' }
#' @export
#'
get_data <- function(table = "LUAD_Academia_protein",
                     action = "expression",
                     genes = c("GAPDH","TNS1")){
  if (length(genes)>1) genes <- paste0(genes,collapse = ",")
  url <- paste0("https://www.jingege.wang/bioinformatics/PCAS/api.php?action=",action,"&table=",table,"&genes=",genes)
  res <- jsonlite::fromJSON(url)
  return(res)
}
