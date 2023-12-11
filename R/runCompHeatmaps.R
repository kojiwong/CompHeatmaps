#' Launch Shiny App for CompHeatmaps
#'
#' A function that launches the Shiny app for CompHeatmaps.
#' The purpose of this app is to allow users to interact and create heatmaps on the fly from raw data.
#' app works. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value, opens up a Shiny page.
#'
#' @examples
#' \dontrun{
#' CompHeatmaps::runCompHeatmaps()
#' }
#'
#' @export
#' @importFrom shiny runApp
#' @import shiny
#' @import colourpicker

runCompHeatmaps <- function() {
  app_dir <- system.file("shiny-scripts",
                        package = "CompHeatmaps")
  action_shiny <- shiny::runApp(app_dir, display.mode = "normal")
  return(action_shiny)
}
# [END]
