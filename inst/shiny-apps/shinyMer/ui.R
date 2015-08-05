# ShinyMer

shiny::fluidPage(
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
          conditionalPanel(condition = "input.simDataType != orig",
                           selectInput("smoothMethod", "Method",
                                       list(names(merMod@frame)))
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
                               value=TRUE),
          shiny::actionButton("goButton", shiny::strong("Run"))
        ),

        shiny::mainPanel(
          shiny::tabsetPanel(type="tabs",
                             shiny::tabPanel("Prediction uncertainty",
                                             shiny::h3("Original function call:"),
                                             # shiny::textOutput("predInt.call"),
                                             shiny::h3("Plot"),
                                             shiny::plotOutput("predIntplot"),
                                             shiny::h3("PredInt data.frame"),
                                             # shiny::downloadButton("downloadData", "Download predict interval data"),
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
    )

