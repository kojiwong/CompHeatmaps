#' Launch Shiny App for CompHeatmaps
#'
#' A function that launches the Shiny app for CompHeatmaps.
#' The purpose of this app is to allow users to interact and create heatmaps on the fly from raw data.
#' app works. The code has been placed in \code{./inst/shiny-scripts}. Uses colourpicker for choosing colours
#' of the low and high values on the heatmap.
#'
#' @return No return value, opens up a Shiny page.
#'
#' @examples
#' \dontrun{
#' CompHeatmaps::runCompHeatmaps()
#' }
#'
#' @references
#' Attali D (2023). colourpicker: A Colour Picker Tool for Shiny and for Selecting Colours in Plots.
#' R package version 1.3.0, <https://CRAN.R-project.org/package=colourpicker>.
#'
#' @export
#' @import shiny
#' @importFrom shiny runApp

runCompHeatmaps <- function() {
  app_dir <- system.file("shiny-scripts",
                        package = "CompHeatmaps")
  action_shiny <- shiny::runApp(app_dir, display.mode = "normal")
  return(action_shiny)
}
# [END]
