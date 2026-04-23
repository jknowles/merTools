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
#' @export
shinyMer <- function(merMod, simData = NULL, pos = 1) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "shinyMer() requires the 'shiny' package. Install it with ",
      "install.packages(\"shiny\").",
      call. = FALSE
    )
  }
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
