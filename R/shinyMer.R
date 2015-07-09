#' Launch a shiny app to explore your merMod interactively
#'
#' \code{shinyMer} launches a shiny app that allows you to interactively
#' explore an estimated merMod using functions from \code{\link{merTools}}.
#'
#' @param merMod An object of class "merMod".
#'
#' @param simData A data.frame to make predictions from (optional). If
#'   NULL, then the user can only make predictions using the data in
#'   the frame slot of the merMod object.
#'
#' @return A shiny app
#'
#' @import shiny
#' @import ggplot2
#' @importFrom DT dataTableOutput
#'
#' @export

shinyMer <- function(merMod, simData=NULL) {
  require(shiny)
  require(ggplot2)

  if (is.null(simData)) {
    df.choices <- c("Model Frame"   = "orig",
                    "Random Obs (NOT IMPLEMENTED YET)"    = "rand",
                    "Average Obs (NOT IMPLEMENTED YET)"   = "mean")
  } else {
    df.choices <- c("User Supplied" = "user",
                    "Model Frame"   = "orig",
                    "Random Obs (NOT IMPLEMENTED YET)"    = "rand",
                    "Average Obs (NOT IMPLEMENTED YET)"   = "mean")

  }

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
                       choices=df.choices,
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
          actionButton("goButton", strong("Run"))
        ),

        mainPanel(
          tabsetPanel(type="tabs",
            tabPanel("Prediction uncertainty",
                     h3("Original function call:"),
                     textOutput("predInt.call"),
                     h3("Plot"),
                     plotOutput("predInt.plot"),
                     h3("PredInt data.frame"),
                     downloadButton("downloadData", "Download predict interval data"),
                     DT::dataTableOutput("predInt.tab")
                     ),
            tabPanel("Parameter Plot",
                     h3("Fixed parameter estimates"),
                     em("Double click on a point to view the estimated values [CI]:"),
                     plotOutput("FEplot", dblclick="FEclick"),
                     h3("Random Parameter estimate"),
                     plotOutput("REplot", click="REclick", width="70%")
                     ),
            tabPanel("Substantive Effects",
                     h3(em("Under construction..."))
                     )
          )
        )
      )
    ),
    #SERVER----
    server = function(input, output, session) {
      rv <- reactiveValues(
          level = NULL,
          scale = NULL,
          n.sims = NULL,
          stat = NULL,
          type = NULL
      )

      sim.df <- eventReactive(input$goButton, {
        if (input$simDataType=="orig") {
          return(data.frame(cbind(merMod@frame, X=1:nrow(merMod@frame))))
        }
        else if (input$simDataType=="user") {
          return(data.frame(cbind(simData, X=1:nrow(simData))))
        }
        else if (input$simDataType=="rand") {
          return("Random Observation is NOT IMPLEMENTED YET")
        }
        else if (input$simDataType=="mean") {
          return("Average Observation is NOT IMPLEMENTED YET")
        }
      })

      observeEvent(input$goButton,
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

      predInt <- reactive({
        predictInterval(merMod,
                        newdata           = sim.df(),
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
          plot.df <- cbind(predInt(), sim.df())

          ytitle <- isolate(ifelse(input$predMetric=="linear.prediction", "Linear Prediction", "Probability"))

          title <- isolate(ifelse(input$simDataType=="user", "User Supplied Data",
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

     output$predInt.tab <- renderDataTable(
       if (input$goButton) {
         data.frame(cbind(sim.df(), predInt()))
       }
     )

     output$downloadData <- downloadHandler(
       filename = "predictIntervalResults.csv",
       content = function(file) {
         write.csv(isolate(cbind(sim.df(), predInt())), file)
       }
     )

     ##Output for Parameter Plot Tab ----
     FEplot.df <- eventReactive(input$goButton, {
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

     clickrows <- reactiveValues(
       clicked = rep(FALSE, length(fixef(merMod)[!grepl("(Intercept)", names(fixef(merMod)), fixed=TRUE)]))
     )

     output$FEplot <- renderPlot(
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

     observeEvent(input$FEclick, {
       clickrows$clicked <- xor(clickrows$clicked, round(input$FEclick$y,1) == as.numeric(FEplot.df()$variable))
     })

     output$REplot <- renderPlot(
       if (input$goButton) {
         scale = qnorm(1-(1-rv$level)/2)
         plotREsim(data=REsim(merMod, rv$n.sims), scale = rv$scale, var=paste(rv$stat, "_eff", sep=""))
       }
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

