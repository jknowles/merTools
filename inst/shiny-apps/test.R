library(merTools)
library(shiny)
library(ggplot2)

merMod <- lmer(Reaction ~ Days + (1|Subject), data = sleepstudy)

shinyApp(ui =
           shinyUI(fluidPage(
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
                 conditionalPanel(condition = "input.simDataType!='orig'",
                                  selectInput("filter", "Filter",
                                              choices = names(merMod@frame))
                                  ),
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
                                      value=TRUE)
               ),
    shiny::mainPanel(
      shiny::tabsetPanel(type="tabs",
              shiny::tabPanel("Prediction uncertainty",
                 textOutput("text1"),
                 dataTableOutput("dt"),
                 plotOutput("predPlot")
               ),
              shiny::tabPanel("Parameters",
                  shiny::h3("Original call"),
                  shiny::h3("Fixed Effects"),
                  plotOutput("feplot")
              ),
              shiny::tabPanel("Substantive Effect",
                  shiny::h3("Effect Sizes")
             )
           )
      )
    )
    )),
          server = shinyServer(function(input, output){
            output$text1 <- renderText({
              paste("You have selected", input$stat)
            })

            predInput <- reactive({
              data <- switch(input$simDataType,
                             "orig" = merMod@frame,
                             "mean" = draw(merMod, type = "average"),
                             "rand" = draw(merMod, type = "random"))
              cbind(predictInterval(merMod, newdata = data,
                              n.sims = input$n.sims, stat = input$stat), data)
            })

            output$dt <- renderDataTable({
              predInput()
            })

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
              FEsim(merMod, n.sims = input$n.sims)
            })

            output$feplot <- renderPlot({
              plotdf <- feData()
              plotdf$y <- ifelse(input$stat == "mean", plotdf$mean_eff,
                                 plotdf$median_eff)
              plotdf$lwr <- plotdf$y - (qnorm(input$alpha/100) * plotdf$sd_eff)
              plotdf$upr <- plotdf$y + (qnorm(input$alpha/100) * plotdf$sd_eff)
              ggplot(plotdf, aes(x = variable, y = y, ymin = lwr, ymax = upr))+
                geom_pointrange() + coord_flip() +
                theme_bw() + theme(axis.text.x = element_blank(),
                                   panel.grid.major.x = element_blank(),
                                   panel.grid.minor.x = element_blank(),
                                   axis.ticks.x = element_blank())

            })
}
)
)
