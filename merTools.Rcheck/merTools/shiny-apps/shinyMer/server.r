#SERVER----
server = function(input, output){
  output$text1 <- renderText({
    paste("You have selected", input$stat)
  })

  predInput <- reactive({
    data <- switch(input$newdataType,
                   "orig" = merMod@frame,
                   "mean" = draw(merMod, type = "average"),
                   "rand" = draw(merMod, type = "random"),
                   "user" = newdata)
    cbind(predictInterval(merMod, newdata = data, level = input$alpha/100,
                          type = input$predMetric,
                          include.resid.var = input$resid.var,
                          n.sims = input$n.sims, stat = input$stat), data)
  })

  output$dt <- renderDataTable({
    predInput()
  })

  output$downloadData <- shiny::downloadHandler(
    filename = "predictIntervalResults.csv",
    content = function(file) {
      write.csv(shiny::isolate(predInput()), file)
    }
  )

  output$predPlot <- renderPlot({
    data <- predInput()
    data$x <- factor(seq(1:nrow(data)))
    ggplot(data, aes(x = x, y = fit, ymin = lwr, ymax = upr)) +
      geom_pointrange() +
      theme_bw() + theme(axis.text.x = element_blank(),
                         panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         axis.ticks.x = element_blank())
  })

  feData <- reactive({
    data <- FEsim(merMod, n.sims = input$n.sims)
    return(data)
  })

  output$feplot <- renderPlot({
    plotdf <- feData()
    scale <- input$alpha/100
    vartmp <- input$stat
    plotFEsim(plotdf, level = scale, stat = vartmp, sd = TRUE,
              intercept = FALSE)
  })

  reData <- reactive({
    data <- REsim(merMod, n.sims = input$n.sims)
    return(data)
  })

  output$replot <- renderPlot({
    plotdf <- reData()
    scale <- input$alpha/100
    vartmp <- input$stat
    plotREsim(plotdf, level = scale, stat = vartmp, sd = TRUE)
  })

  output$call <- renderPrint({
    merMod@call
  })

  reEffInput <- reactive({
    data <- switch(input$newdataType,
                   "orig" = merMod@frame,
                   "mean" = draw(merMod, type = "average"),
                   "rand" = draw(merMod, type = "random"),
                   "user" = newdata)
    if(nrow(data) > 12){
      warning("Too much data selected, only using top 12 rows.")
      data <- data[1:12, ]
    }
    return(data)
  })

  groupData <- reactive({
    plotdf <- REimpact(merMod, newdata = reEffInput(),
                       groupFctr = input$group,
                       term = input$term,
                       level = input$alpha/100,
                       breaks = input$nbin,
                       type = input$predMetric,
                       include.resid.var = input$resid.var,
                       n.sims = input$n.sims, stat = input$stat)
    plotdf$upr <- qnorm(input$alpha/100) * plotdf$AvgFitSE
    plotdf$lwr <- qnorm(input$alpha/100) * plotdf$AvgFitSE
    plotdf$upr <- plotdf$AvgFit + plotdf$upr
    plotdf$lwr <- plotdf$AvgFit - plotdf$lwr
    plotdf$bin <- factor(plotdf$bin)
    return(plotdf)
  })

  output$gPlot <- renderPlot({
    ggplot(groupData(), aes(x = bin, y = AvgFit, ymin = lwr, ymax = upr)) +
      geom_pointrange() + facet_wrap(~case) +
      theme_bw() + labs(x = "Bin", y = "Value of DV",
                        title = "Impact of grouping term for selected case")
  })

  wiggleData <- reactive({
    valLookup <- unique(merMod@frame[, input$fixef])
    if(class(valLookup) %in% c("numeric", "integer")){
      newvals <- seq(min(valLookup), max(valLookup), length.out = 20)
    } else{
      if(length(valLookup) < 50){
        newvals <- valLookup
      } else{
        newvals <- sample(valLookup, 50)
      }
    }
    plotdf <- wiggle(reEffInput(), input$fixef, values = newvals)
    plotdf <- cbind(plotdf, predictInterval(merMod, newdata=plotdf,
                                            type = input$predMetric,
                                            level = input$alpha/100,
                                            include.resid.var = input$resid.var,
                                            n.sims = input$n.sims, stat = input$stat))
    plotdf$X <- plotdf[, input$fixef]
    plotdf$case <- rep(1:length(newvals), length = nrow(reEffInput()))
    return(plotdf)
  })

  output$re.ui <- renderUI({
    choices <- names(ranef(merMod)[[input$group]])
    selectInput("term", "Group Term:",
                choices = choices,
                selected = choices[1])
  })

  output$wigglePlot <- renderPlot({
    ggplot(wiggleData(), aes(x = X, y = fit, ymin = lwr,
                             ymax = upr)) +
      geom_pointrange() + facet_wrap(~case) +
      theme_bw() + labs(y = "Simulated Value of DV",
                        title = "Impact of selected fixed effect for
                        selected cases.")
  })

  }

