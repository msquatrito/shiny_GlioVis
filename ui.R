library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)

shinyUI(fluidPage(
  # progressInit() must be called somewhere in the UI in order for the progress UI to actually appear
  progressInit(),

  sidebarLayout(
                
                sidebarPanel(
                  #       tags$head(tags$style(".well {background-color: black; }")),
                  img(src = "GlioVis_logo.jpg", height = 90, width = 270),
                  br(),
                  br(),
                  selectInput("dataset", 
                              label = h4("Dataset"),
                              choices = c("TCGA GBM", "TCGA Lgg", "Rembrandt", 
                                          "Gravendeel", "Phillips", "Murat", "Freije"),
                              selected = "TCGA GBM"),
                  br(),
                  selectInput("plotTypeSel", 
                              label = h4("Plot type"), choices = ""),
                  br(),
                  selectizeInput("gene", 
                                 label = h4("Gene"), choices = NULL, selected = NULL,
                                 options = list(placeholder = "Enter gene, eg: EGFR", 
                                                plugins = list('restore_on_backspace'))),
                  br(),
                  uiOutput("help")
                ),
                
                
                mainPanel(
                  
                  h3(textOutput("caption")),
                  
                  tabsetPanel(
                    
                    tabPanel("Plot", icon = icon("bar-chart-o"),
                             
                             fluidRow(
                               column(4, 
                                      wellPanel(
                                        strong("Plot options:"),
                                        checkboxInput("scale", "Scale y axis", FALSE),
                                        checkboxInput("colStrip", "Color stripchart", FALSE),
                                        checkboxInput("colBox", "Color box", FALSE))),
                               
                               column(4, 
                                      wellPanel(
                                        strong("Statistic:"),
                                        checkboxInput("stat", "Tukey's HSD (plot)", FALSE),
                                        checkboxInput("statTable", "Tukey's HSD (table)", FALSE)))),
                             
                             plotOutput("plot"),
                             br(),
                             conditionalPanel(
                               condition = "input.statTable",
                               strong(" Tukey's Honest Significant Difference (HSD):"),
                               verbatimTextOutput("tukeyTest")),
                             br(),
                             
                             wellPanel( 
                               selectInput(
                                 inputId = "downloadPlotFileType",
                                 label   = h5("Select download file type"),
                                 choices = list(
                                   "PDF"  = "pdf",
                                   "BMP"  = "bmp",
                                   "JPEG" = "jpeg",
                                   "PNG"  = "png")),
                               
                               # Allow the user to set the height and width of the plot download.
                               h5(HTML("Set download image dimensions<br>(units are inches for PDF, pixels for all other formats)")),
                               
                               numericInput(
                                 inputId = "downloadPlotHeight",
                                 label = "Height (inches)",
                                 value = 7,
                                 min = 1,
                                 max = 100),
                               
                               numericInput(
                                 inputId = "downloadPlotWidth",
                                 label = "Width (inches)",
                                 value = 7,
                                 min = 1,
                                 max = 100)
                             ),
                             br(),               
                             downloadButton('downloadPlot', 'Download')),
                    
                    tabPanel("Survival", icon = icon("user-md"),
                             tabsetPanel(
                               tabPanel("Km plot",
                                        fluidRow(
                                          column(4, 
                                                 selectInput("histologySurv", strong("Histology:"), 
                                                             c("GBM"))
                                          ),
                                          column(4,
                                                 selectInput("cutoff", strong("Select cutoff:"),
                                                             c("median", "lower quartile", "upper quartile", "quartiles"))
                                          ),
                                          column(4, 
                                                 selectInput("subtypeSurv", strong("Subtype (GBM only):"), 
                                                             c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                                          )
                                        ),
                                        conditionalPanel(
                                          condition = "input.histologySurv == 'GBM'",
                                          checkboxInput("gcimpSurv", "Exclude G-CIMP samples", FALSE)
                                        ),
                                        plotOutput("survPlot"),
                                        downloadButton('downloadsurvPlot', 'Download')
                               ),
                               
                               tabPanel("HR plot", 
                                        wellPanel(
                                          helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                                          helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                              to update the survival plot. The blue line represents the current selection.")),
                                          br(),
                                          checkboxInput("quantile", "Show quantiles", TRUE),
                                          checkboxInput("gcimp", "Exclude G-CIMP samples", FALSE),
                                          plotOutput("hazardPlot", clickId="densityClick")
                                        ),
                                        plotOutput("kmPlot"),
                                        downloadButton('downloadkmPlot', 'Download')
                               )
                             )
                    ),
                    
                    
                    tabPanel("Correlations", icon = icon("list-alt"),
                             tabsetPanel(
                               tabPanel("Correlation plot", # icon = icon("bar-chart-o"),
                                        fluidRow(
                                          column(4,
                                                 selectInput("gene1", h5("Gene 1"), "", selectize = TRUE)
                                          ),
                                          column(4, 
                                                 selectInput("gene2", h5("Gene 2"), "", selectize = TRUE)
                                          )),
                                        br(),
                                        fluidRow(
                                          column(3, 
                                                 selectInput("histologyCorr", h5("Histology:"), choices = "")
                                          ),
                                          column(3, 
                                                 selectInput("subtype", h5("Subtype (GBM):"), 
                                                             c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                                          ),
                                          column(3, inputPanel(
                                            radioButtons("colorBy", 
                                                         label = strong("Color by:"), 
                                                         choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"))
                                          )),
                                          column(3, inputPanel(
                                            radioButtons("separateBy", 
                                                         label = strong("Separate by:"), 
                                                         choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"))
                                          ))),
                                        br(),
                                        plotOutput("corrPlot"),
                                        downloadButton('downloadcorrPlot', 'Download plot'), 
                                        br(),
                                        verbatimTextOutput("corrTest"), id = "corrPlot") ,
                               
                               tabPanel("Correlation data", # icon = icon("list-alt"),
                                        fluidRow(
                                          column(4, 
                                                 selectInput("histologyCorrTable", "Histology:", choices = "")
                                          ),
                                          column(4, 
                                                 selectInput("cor", "Correlation:", 
                                                             c("All", "Positive", "Negative"))),
                                          column(4, 
                                                 selectInput("sign", "Signficance:", 
                                                             c(0.05, 0.01)))),
                                        br(),
                                        downloadButton('downloadCorrData', 'Download data'),
                                        dataTableOutput("corrData"), id = "corrData")
                             )
                    ),
                    
                    tabPanel("Data", icon = icon("table"),
                             uiOutput("piePlots"),
                             fluidRow(
                               column(12,
                                      verbatimTextOutput("summary"))),
                             br(),
                             downloadButton('downloadData', 'Download table'), 
                             br(),
                             dataTableOutput("table")),
                    
                    tabPanel("About", includeMarkdown("tools/about.Rmd"))
                    
                  )
                )
  )
))
