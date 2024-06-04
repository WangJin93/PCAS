#' Run PCAS Shiny App
#' @import shiny
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' PCAS_app()
#' }
PCAS_app <- function() {
  shiny::shinyAppFile(system.file("shinyapp", "app.R", package = "PCAS"))
}
