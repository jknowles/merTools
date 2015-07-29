#SERVER----

library(merTools)
library(lme4)
library(ggplot2)

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

  shiny::observeEvent(input$goButton > 0,
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


  subplot.df <- shiny::eventReactive(input$goButton > 0, {
    tmp <- sim.df()
    subplot.df <- predictInterval(merMod, newdata = tmp,
                                 level = rv$level, n.sims = rv$n.sims,
                                 stat = rv$stat, type = rv$type,
                                 include.resid.var = rv$include.resid.var)
    outdata <- cbind(tmp, subplot.df)
    return(outdata)
  }
  )

  output$predInt.call <- shiny::renderPrint(
    getCall(merMod)
  )

  output$predIntplot <- shiny::renderPlot({
#       title <- shiny::isolate(ifelse(input$predMetric=="linear.prediction", "Linear Prediction", "Probability"))
#
#       ytitle <- shiny::isolate(ifelse(input$simDataType=="user", "User Supplied Data",
#                                      ifelse(input$simDataType=="orig", "Model Frame Data",
#                                             ifelse(input$simDataType=="rand", "Random Observation(s)",
#                                                    ifelse(input$simDataType=="mean", "Averge Observation(s)")))))
      plot <-  ggplot(aes(x=X, y=fit, ymin=lwr, ymax=upr), data=   subplot.df()) +
        geom_errorbar(color="gray50") +
        geom_point() +
        theme_bw() +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks.x=element_blank(), axis.line.x=element_blank()) # +
#         labs(y=ytitle, title=paste("Predicted values for", title))
      #Highlight selected rows in plot
#       s <- input$predInt.tab_rows_selected
#       if (length(s)) {
#         return(plot +
#           geom_errorbar(data=plot.df[s,], color="red") +
#           geom_point(data=plot.df[s,], color="red"))
#       } else {
#         return(plot)
#       }
      return(plot)

  })

  output$predInt.tab <- shiny::renderDataTable({
    return(subplot.df())
  })

  output$downloadData <- shiny::downloadHandler(
    filename = "predictIntervalResults.csv",
    content = function(file) {
      write.csv(shiny::isolate(cbind(sim.df(), predInt())), file)
    }
  )

  ##Output for Parameter Plot Tab ----
  FEplot.df <- shiny::eventReactive(input$goButton > 0, {
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
    if (input$goButton > 0) {
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
    if (input$goButton > 0) {
      scale = qnorm(1-(1-rv$level)/2)
      plotREsim(data=REsim(merMod, rv$n.sims), scale = rv$scale, var=paste(rv$stat, "_eff", sep=""))
    }
  )
}
