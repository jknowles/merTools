#SERVER----
server = function(input, output, session){
  output$text1 <- renderText({
    paste("You have selected", input$stat)
  })

  # Optional subset of the model frame for the Random/Average draws (issue #32).
  # Returns a named list suitable for draw(..., varList = ) or NULL when the user
  # has not chosen a variable to restrict on.
  subsetVarList <- reactive({
    if (is.null(input$subsetVar) || identical(input$subsetVar, "") ||
        is.null(input$subsetVal) || identical(input$subsetVal, "")) {
      return(NULL)
    }
    stats::setNames(list(input$subsetVal), input$subsetVar)
  })

  # Value picker for the chosen subset variable, populated from the model frame.
  output$subsetVal_ui <- renderUI({
    if (is.null(input$subsetVar) || identical(input$subsetVar, "")) {
      return(NULL)
    }
    vals <- sort(unique(as.character(merMod@frame[, input$subsetVar])))
    selectInput("subsetVal",
                paste0("Value of ", input$subsetVar, ":"),
                choices = vals,
                selected = vals[1])
  })

  # Resolve the requested scenario into a data.frame of cases to predict from.
  drawData <- function() {
    switch(input$newdataType,
           "orig" = merMod@frame,
           "mean" = draw(merMod, type = "average", varList = subsetVarList()),
           "rand" = draw(merMod, type = "random", varList = subsetVarList()),
           "user" = newdata)
  }

  predInput <- reactive({
    data <- drawData()
    cbind(predictInterval(merMod, newdata = data, level = input$alpha/100,
                          type = input$predMetric,
                          include.resid.var = input$resid.var,
                          n.sims = input$n.sims, stat = input$stat), data)
  })

  if ("DT" %in% rownames(installed.packages())) {
    output$dt <- DT::renderDT({
      predInput()
    })
  } else {
    output$dt <- renderTable({
      predInput()
    })
  }

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

  # ---- Model Summary tab (issue #78) ----
  output$modelInfo <- renderTable({
    merTools::modelInfo(merMod)
  }, digits = 3)

  output$summaryCall <- renderPrint({
    print(merMod@call)
  })

  output$fixefTab <- renderTable({
    fe <- lme4::fixef(merMod)
    data.frame(term = names(fe), estimate = as.numeric(fe),
               check.names = FALSE)
  }, digits = 4)

  output$vcTab <- renderPrint({
    print(lme4::VarCorr(merMod), comp = c("Variance", "Std.Dev."))
  })

  output$ngrpsTab <- renderTable({
    ng <- lme4::ngrps(merMod)
    data.frame(grouping_factor = names(ng), n_levels = as.integer(ng),
               check.names = FALSE)
  })

  reEffInput <- reactive({
    data <- drawData()
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
    # Two-sided interval half-width from the requested level.
    z <- qnorm(1 - (1 - input$alpha/100) / 2)
    plotdf$upr <- plotdf$AvgFit + z * plotdf$AvgFitSE
    plotdf$lwr <- plotdf$AvgFit - z * plotdf$AvgFitSE
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
    base <- reEffInput()
    valLookup <- unique(merMod@frame[, input$fixef])
    if (is.numeric(valLookup)) {
      newvals <- seq(min(valLookup), max(valLookup), length.out = 20)
    } else if (length(valLookup) < 50) {
      newvals <- valLookup
    } else {
      newvals <- sample(valLookup, 50)
    }
    plotdf <- wiggle(base, input$fixef, values = list(newvals))
    plotdf <- cbind(plotdf, predictInterval(merMod, newdata = plotdf,
                                            type = input$predMetric,
                                            level = input$alpha/100,
                                            include.resid.var = input$resid.var,
                                            n.sims = input$n.sims, stat = input$stat))
    plotdf$X <- plotdf[, input$fixef]
    # wiggle() stacks the base data once per value, so each base observation
    # (the "case" we facet on) repeats across the value blocks.
    plotdf$case <- rep(seq_len(nrow(base)), times = length(newvals))
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
