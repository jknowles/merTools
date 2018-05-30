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
#' @param pos The position of the environment to export function arguments to.
#' Defaults to 1, the global environment, to allow shiny to run.
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
#' @importFrom shiny runApp
#' @importFrom DT dataTableOutput
#' @export

shinyMer <- function(merMod, simData = NULL, pos = 1) {
  envir = as.environment(pos)
  if(exists("simData")){
    expParm <- function(x, y) assign(".shinyMerPar", list("merMod" = x, "simData" = y), envir = envir)
    expParm(x = merMod, y = simData)
  } else{
    expParm2 <- function(x) assign(".shinyMerPar", list("merMod" = x, "simData" = NULL), envir = envir)
    expParm2(x = merMod)
  }
  appDir <- system.file("shiny-apps", "shinyMer", package = "merTools")
  shiny::runApp(appDir, display.mode = "normal")
}
