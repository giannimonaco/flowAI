library(shiny)

shinyUI(fluidPage(
    titlePanel("Interactive Quality Control for Flow Cytometry Data"),

    fluidRow(
        column(3,

               h4("Input:"),
               fileInput('fcsFiles', strong('Choose FCS file:'), multiple = FALSE,
                         accept = c('text/fcs', '.fcs')),
               actionButton("goButton", "Submit!"),

               hr(),
               h4("Summary:"),
               textOutput("summaryText1"),
               textOutput("summaryText2"),

               hr(),
               h4("Parameters:"),
               numericInput("timeLenth", label = h5("Time step (sec)"), value = 0.1, step = 0.1),
               uiOutput("signalBinSize"),

               hr(),
               h4("Output:"),
               downloadButton('downloadGoodFCS', 'Download High Q file (FCS)'),
               br(),
               downloadButton('downloadBadFCS', 'Download Low Q file (FCS)'),

               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ",
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
               tabsetPanel(type = "pills",

                           tabPanel("Flow Rate", fluidPage(
                               hr(),
                               textOutput("flowRateSummary"),
                               hr(),
                               plotOutput("flowRatePlot"),
                               hr(),
                               fluidRow(
                                   column(4, offset = 1,
                                          uiOutput("timeSlider")
                                   ),
                                   column(4, offset = 2,
                                          uiOutput("rateSlider")
                                   )
                               )
                           )),

                           tabPanel("Signal Acquisition", fluidPage(
                               hr(),
                               textOutput("flowSignalSummary"),
                               hr(),
                               uiOutput("signalBinSlider"),
                               hr(),
                               plotOutput("flowSignalPlot", height = "800px")
                           )),

                           tabPanel("Dynamic Range", fluidPage(
                               hr(),
                               fluidRow(
                                   column(5,
                                          textOutput("flowMarginSummary")
                                   ),
                                   column(3, offset = 2,
                                          checkboxInput("checkbox", label = "Apply Margins Check", value = TRUE)
                                   )
                               ),
                               hr(),
                               plotOutput("flowMarginPlot")
                            ))
               )
        )
    )
))
