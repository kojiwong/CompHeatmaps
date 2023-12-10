#' Launch Shiny App for CompHeatmaps
#'
#' A function that launches the Shiny app for CompHeatmaps.
#' The purpose of this app is to allow users to interact and create heatmaps on the fly from raw data.
#' app works. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' CompHeatmaps::runCompHeatmaps()
#' }
#'
#' @export
#' @importFrom shiny runApp

runCompHeatmaps <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "CompHeatmaps")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
