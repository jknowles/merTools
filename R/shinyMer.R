utils::globalVariables(c("X", "fit", "lwr", "upr", "variable", "lci", "uci", "label"))
#' Launch a shiny app to explore your merMod interactively
#'
#' \code{shinyMer} launches a shiny app that allows you to interactively
#' explore an estimated merMod using functions from \code{merTools}.
#'
#' @param merMod An object of class "merMod".
#'
#' @param simData A data.frame to make predictions from (optional). If
#'   NULL, then the user can only make predictions using the data in
#'   the frame slot of the merMod object.
#'
#' @return A shiny app
#'
#' @import ggplot2
#' @importFrom shiny shinyApp
#' @importFrom shiny fluidPage
#' @importFrom shiny titlePanel
#' @importFrom shiny sidebarLayout
#' @importFrom shiny sidebarPanel
#' @importFrom shiny radioButtons
#' @importFrom shiny numericInput
#' @importFrom shiny checkboxInput
#' @importFrom shiny actionButton
#' @importFrom shiny mainPanel
#' @importFrom shiny tabsetPanel
#' @importFrom shiny tabPanel
#' @importFrom shiny h3
#' @importFrom shiny textOutput
#' @importFrom shiny plotOutput
#' @importFrom shiny downloadButton
#' @importFrom shiny em
#' @importFrom shiny reactiveValues
#' @importFrom shiny eventReactive
#' @importFrom shiny observeEvent
#' @importFrom shiny reactive
#' @importFrom shiny renderPrint
#' @importFrom shiny renderPlot
#' @importFrom shiny isolate
#' @importFrom shiny renderPrint
#' @importFrom shiny downloadHandler
#' @importFrom shiny strong
#' @importFrom DT dataTableOutput
#' @export
runExample <- function(merMod, simData = NULL) {
  if (is.null(simData)) {
    df.choices <- c("Model Frame"   = "orig",
                    "Random Obs"    = "rand",
                    "Average Obs"   = "mean")
  } else {
    df.choices <- c("User Supplied" = "user",
                    "Model Frame"   = "orig",
                    "Random Obs"    = "rand",
                    "Average Obs"   = "mean")

  }
  appDir <- system.file("shiny-apps", "shinyMer", package = "mypackage")
  shiny::runApp(appDir, display.mode = "normal")
}
