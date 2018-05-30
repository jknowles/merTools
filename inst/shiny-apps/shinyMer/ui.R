# ShinyMer
shinyUI(fluidPage(
  shiny::titlePanel("Explore your merMod interactively"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::radioButtons("newdataType",
                          "Simulated data scenario",
                          choices=df.choices,
                          selected=NULL),
      #                  conditionalPanel(condition = "input.newdataType!='orig'",
      #                                   selectInput("filter", "Filter",
      #                                               choices = names(merMod@frame))
      #                                   ),
      shiny::numericInput("n.sims",
                          label="Simulations (Max=1,000)",
                          value=100,
                          min=1,
                          max=1000),
      conditionalPanel(condition="input.conditionedPanels==3",
                       helpText("Here you can compare impact of
                                changing input variables on the
                                outcome variable for selected cases."),
                       selectInput("group", "Group Factor:",
                                   choices = names(ranef(merMod)),
                                   selected = NULL),
                       uiOutput("re.ui"),
                       sliderInput("nbin", "Effect Bins", min = 3,
                                   max = 10, value = 4, step = 1),
                       helpText("And modify your fixed effects"),
                       selectInput("fixef", "Fixed Effect:",
                                   choices = all.vars(nobars(formula(merMod)))[-1],
                                   selected = NULL)
                       ),
      shiny::numericInput("alpha",
                          label="Credible Interval (%)",
                          value=95,
                          min=0,
                          max=100),
      shiny::radioButtons("stat",
                          "Measure of central tendency",
                          choices=c("Median"="median", "Mean"="mean"),
                          selected=NULL),
      shiny::radioButtons("predMetric",
                          "Prediction metric",
                          choices=c("Linear Predictor"="linear.prediction",
                                    "Probability"="probability"),
                          selected=NULL),
      shiny::checkboxInput("resid.var",
                           label="Include Residual Variation",
                           value=TRUE)
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(type="tabs",
                         shiny::tabPanel("Prediction uncertainty",
                                         shiny::h3("Prediction Intervals:"),
                                         plotOutput("predPlot"),
                                         shiny::h3("All Predictions"),
                                         dataTableOutput("dt"),
                                         shiny::downloadButton("downloadData", "Download predict interval data")
                                         ,value = 1
                         ),
                         shiny::tabPanel("Parameters",
                                         shiny::h3("Original call"),
                                         textOutput("call"),
                                         shiny::h3("Fixed Effects"),
                                         plotOutput("feplot"),
                                         shiny::h3("Group Effects"),
                                         plotOutput("replot"), value = 2
                         ),
                         shiny::tabPanel("Substantive Effect",
                                         shiny::h3("Effect Sizes"),
                                         plotOutput("gPlot"),
                                         shiny::h3("Fixef Effect Impact"),
                                         plotOutput("wigglePlot"), value = 3
                         ), id = "conditionedPanels"
      )
    )
)
  ))
