#' Launch a shiny app to explore your merMod interactively
#'
#' \code{shinyMer} launches a shiny app that allows you to interactively
#' explore an estimated merMod using functions from \code{\link{merTools}}.
#'
#' @param merMod An object of class "merMod".
#'
#' @return A shiny app
#'
#' @import shiny
#' @import ggplot2
#'
#' @export

shinyMer <- function(merMod, simData=NULL) {
  require(shiny)
  require(ggplot2)

  shinyApp(
    #UI----
    ui = fluidPage(
      titlePanel("Explore your merMod interactively"),

      sidebarLayout(
        sidebarPanel(
          radioButtons("stat",
                       "Measure of central tendency",
                       choices=c("Median"="median", "Mean"="mean"),
                       selected=NULL),
          radioButtons("predMetric",
                       "Prediction metric",
                       choices=c("Linear Predictor"="linear.prediction",
                                 "Probability"="probability"),
                       selected=NULL),
          radioButtons("simDataType",
                       "Simulated data scenario",
                       choices=c("Model Frame"   = "orig",
                                 "User Supplied" = "user",
                                 "Random Obs"    = "rand",
                                 "Average Obs"   = "mean"),
                       selected=NULL),
          numericInput("n.sims",
                       label="Simulations (Max=10,000)",
                       value=20,
                       min=1,
                       max=10000),
          numericInput("alpha",
                       label="Credible Interval (%)",
                       value=95,
                       min=0,
                       max=100),
          checkboxInput("resid.var",
                        label="Include Residual Variation",
                        value=TRUE),
          actionButton("goButton", "Go!")
        ),

        mainPanel(
          tabsetPanel(type="tabs",
            tabPanel("Random Effects Uncertainty",
                     em("Original function call:"),
                     textOutput("predInt.call"),
                     plotOutput("predInt")
                     )
          )
        )
      )
    ),
    #SERVER----
    server = function(input, output, session) {
      output$predInt.call <- renderPrint(
        getCall(merMod)
      )
      output$predInt <- renderPlot({
        input$goButton
        #check simData type
        if (input$simDataType=="orig") {
          simData=merMod@frame
        } else if (input$simDataType=="user") {
          if (is.null(simData)) {
            error("You failed to supply simData")
          }
          simData=simData
        } else if (input$simDataType=="rand") {
          error("Not Yet Implemented")
        } else if (input$simDataType=="mean") {
          error("Not Yet Implemented")
        }
        isolate({
          predInt <- predictInterval(merMod,
                                     newdata=simData,
                                     level=input$alpha/100,
                                     n.sims=input$n.sims,
                                     stat=input$stat,
                                     type=input$predMetric,
                                     include.resid.var=input$resid.var)

          ggplot(x=1:nrow(predInt), y=fit, ymin=lci, ymax=uci, data=predInt) +
            geom_errorbar(color="grey50") + geom_point()
        })
      })
    }
  )
}


#scraps
##Extract levels from model for input buttons
#lvls <- NULL
#for (i in 1:getME(merMod, "n_rfacs")) {
#  lvl.name <- names(ranef(merMod))[i]
#  lvls <- c(lvls, paste(lvl.name, ": ", ifelse(names(ranef(merMod)[[i]])=="(Intercept)", "Intercept", names(ranef(merMod)[[i]])), sep=""))
#}
#radioButtons("level",
#             "Choose random term to display",
#             choices=lvls,
#             selected=NULL),

