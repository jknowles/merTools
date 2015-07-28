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
shinyMer <- function(merMod, simData=NULL) {
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

  shiny::shinyApp(
    #UI----
    ui = shiny::fluidPage(
      shiny::titlePanel("Explore your merMod interactively"),

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::radioButtons("stat",
                       "Measure of central tendency",
                       choices=c("Median"="median", "Mean"="mean"),
                       selected=NULL),
          shiny::radioButtons("predMetric",
                       "Prediction metric",
                       choices=c("Linear Predictor"="linear.prediction",
                                 "Probability"="probability"),
                       selected=NULL),
          shiny::radioButtons("simDataType",
                       "Simulated data scenario",
                       choices=df.choices,
                       selected=NULL),
          shiny::numericInput("n.sims",
                       label="Simulations (Max=10,000)",
                       value=100,
                       min=1,
                       max=10000),
          shiny::numericInput("alpha",
                       label="Credible Interval (%)",
                       value=95,
                       min=0,
                       max=100),
          shiny::checkboxInput("resid.var",
                        label="Include Residual Variation",
                        value=TRUE),
          shiny::actionButton("goButton", shiny::strong("Run"))
        ),

        shiny::mainPanel(
          shiny::tabsetPanel(type="tabs",
            shiny::tabPanel("Prediction uncertainty",
                     shiny::h3("Original function call:"),
                     shiny::textOutput("predInt.call"),
                     shiny::h3("Plot"),
                     shiny::plotOutput("predInt.plot"),
                     shiny::h3("PredInt data.frame"),
                     shiny::downloadButton("downloadData", "Download predict interval data"),
                     DT::dataTableOutput("predInt.tab")
                     ),
            shiny::tabPanel("Parameter Plot",
                     shiny::h3("Fixed parameter estimates"),
                     shiny::em("Double click on a point to view the estimated values [CI]:"),
                     shiny::plotOutput("FEplot", dblclick="FEclick"),
                     shiny::h3("Random Parameter estimate"),
                     shiny::plotOutput("REplot", click="REclick", width="70%")
                     ),
            shiny::tabPanel("Substantive Effects",
                     shiny::h3(shiny::em("Under construction..."))
                     )
          )
        )
      )
    ),
    #SERVER----
    server = function(input, output, session) {
      rv <- shiny::reactiveValues(
          level = NULL,
          scale = NULL,
          n.sims = NULL,
          stat = NULL,
          type = NULL
      )

      sim.df <- shiny::eventReactive(input$goButton, {
        if (input$simDataType=="orig") {
          return(data.frame(cbind(merMod@frame, X=1:nrow(merMod@frame))))
        }
        else if (input$simDataType=="user") {
          return(data.frame(cbind(simData, X=1:nrow(simData))))
        }
        else if (input$simDataType=="rand") {
          return(data.frame(cbind(draw(merMod, type = "random"), X = 1)))
        }
        else if (input$simDataType=="mean") {
          return(data.frame(cbind(draw(merMod, type = "average"), X =1 )))
        }
      })

      shiny::observeEvent(input$goButton,
                    {
                      rv$level  = input$alpha/100
                      rv$scale  = qnorm(1-(1-(input$alpha/100))/2)
                      rv$n.sims = input$n.sims
                      rv$stat   = input$stat
                      rv$type   = input$predMetric
                      rv$include.resid.var = input$resid.var
                    }
                  )

      ##Output for Uncertainty tab ----

      predInt <- shiny::reactive({
        predictInterval(merMod,
                        newdata           = sim.df(),
                        level             = rv$level,
                        n.sims            = rv$n.sims,
                        stat              = rv$stat,
                        type              = rv$type,
                        include.resid.var = rv$include.resid.var)
      })

      output$predInt.call <- shiny::renderPrint(
        getCall(merMod)
      )

      output$predInt.plot <- shiny::renderPlot({
        if (input$goButton) {
          plot.df <- cbind(predInt(), sim.df())

          ytitle <- shiny::isolate(ifelse(input$predMetric=="linear.prediction", "Linear Prediction", "Probability"))

          title <- shiny::isolate(ifelse(input$simDataType=="user", "User Supplied Data",
                           ifelse(input$simDataType=="orig", "Model Frame Data",
                           ifelse(input$simDataType=="rand", "Random Observation(s)",
                           ifelse(input$simDataType=="mean", "Averge Observation(s)")))))

         plot <-  ggplot(aes(x=X, y=fit, ymin=lwr, ymax=upr), data=plot.df) +
                    geom_errorbar(color="gray50") +
                    geom_point() +
                    theme_bw() +
                    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(), axis.line.x=element_blank()) +
                    labs(y=ytitle, title=paste("Predicted values for", title))
          #Highlight selected rows in plot
          s <- input$predInt.tab_rows_selected
          if (length(s)) {
            plot +
              geom_errorbar(data=plot.df[s,], color="red") +
              geom_point(data=plot.df[s,], color="red")
          } else {
            plot
          }

        }
     })

     output$predInt.tab <- shiny::renderDataTable(
       if (input$goButton) {
         data.frame(cbind(sim.df(), predInt()))
       }
     )

     output$downloadData <- shiny::downloadHandler(
       filename = "predictIntervalResults.csv",
       content = function(file) {
         write.csv(shiny::isolate(cbind(sim.df(), predInt())), file)
       }
     )

     ##Output for Parameter Plot Tab ----
     FEplot.df <- shiny::eventReactive(input$goButton, {
       FEplot.df <- FEsim(merMod, rv$n.sims)
       FEplot.df$fit <- FEplot.df[,paste(rv$stat, "_eff", sep="")]
       FEplot.df$uci <- FEplot.df$fit + rv$scale*FEplot.df$sd_eff
       FEplot.df$lci <- FEplot.df$fit - rv$scale*FEplot.df$sd_eff
       FEplot.df$label <- paste(
         formatC(FEplot.df$fit, digits=2), " [",
         formatC(FEplot.df$lci, digits=2), ", ",
         formatC(FEplot.df$uci, digits=2), "]", sep=""
       )
       FEplot.df <- subset(FEplot.df, variable!="(Intercept)")
       FEplot.df$variable <- factor(FEplot.df$variable)
       return(FEplot.df)
     }
     )

     clickrows <- shiny::reactiveValues(
       clicked = rep(FALSE, length(fixef(merMod)[!grepl("(Intercept)", names(fixef(merMod)), fixed=TRUE)]))
     )

     output$FEplot <- shiny::renderPlot(
       if (input$goButton) {
         FEplot <- ggplot(aes(x=variable, y=fit, ymin=lci, ymax=uci, label=label), data=FEplot.df()) +
                     geom_errorbar(width=0.2) +
                     geom_point(size=I(3)) +
                     geom_hline(yintercept = 0, color = "red") +
                     theme_bw() + coord_flip()
         if (sum(clickrows$clicked)!=0) {
           return(FEplot + geom_text(data=FEplot.df()[clickrows$clicked,], vjust=1))
         } else {
           return(FEplot)
         }
       }
     )

     shiny::observeEvent(input$FEclick, {
       clickrows$clicked <- xor(clickrows$clicked, round(input$FEclick$y,1) == as.numeric(FEplot.df()$variable))
     })

     output$REplot <- shiny::renderPlot(
       if (input$goButton) {
         scale = qnorm(1-(1-rv$level)/2)
         plotREsim(data=REsim(merMod, rv$n.sims), scale = rv$scale, var=paste(rv$stat, "_eff", sep=""))
       }
     )
   }
  )
}

