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
                       choices=c("User Supplied" = "user",
                                 "Model Frame"   = "orig",
                                 "Random Obs (NOT IMPLEMENTED YET)"    = "rand",
                                 "Average Obs (NOT IMPLEMENTED YET)"   = "mean"),
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
            tabPanel("Prediction uncertainty",
                     h3("Original function call:"),
                     textOutput("predInt.call"),
                     h3("Plot"),
                     plotOutput("predInt.plot"),
                     h3("PredInt data.frame"),
                     dataTableOutput("predint.tab")
                     )
          )
        )
      )
    ),
    #SERVER----
    server = function(input, output, session) {
      rv <- reactiveValues(
          level = NULL,
          n.sims = NULL,
          stat = NULL,
          type = NULL
      )

      observeEvent(input$goButton,
                    {
                      rv$level  = input$alpha/100
                      rv$n.sims = input$n.sims
                      rv$stat   = input$stat
                      rv$type   = input$predMetric
                      rv$include.resid.var = input$resid.var
                      if (input$simDataType=="orig") {
                        simData <- merMod@frame
                      }
                      else if (input$simDataType=="user") {
                        simData <- simData
                      }
                      else if (input$simDataType=="rand") {
                        simData <- "Random Observation is NOT IMPLEMENTED YET"
                      }
                      else if (input$simDataType=="rand") {
                        simData <- "Average Observation is NOT IMPLEMENTED YET"
                      }
                    }
                  )

      predInt <- reactive({
        predictInterval(merMod,
                        newdata           = simData,
                        level             = rv$level,
                        n.sims            = rv$n.sims,
                        stat              = rv$stat,
                        type              = rv$type,
                        include.resid.var = rv$include.resid.var)
      })

      output$predInt.call <- renderPrint(
        getCall(merMod)
      )

      output$predInt.plot <- renderPlot({
        if (input$goButton) {
          predInt <- predInt()

          predInt$x <- 1:nrow(predInt)

          ytitle <- isolate(ifelse(input$predMetric=="linear.prediction", "Linear Prediction", "Probability"))

          ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr), data=predInt) +
            geom_errorbar(color="gray50") +
            geom_point() +
            theme_bw() +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(), axis.line.x=element_blank()) +
            labs(y=ytitle)
        }
     })

     output$predInt.tab <- renderDataTable(
       as.data.frame(cbind(simData, predInt()))
     )

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

