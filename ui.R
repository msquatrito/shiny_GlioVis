library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)
library (shinysky)
library (shinyBS)

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
                  selectizeInput("gene", 
                                 label = h4("Gene"), choices = NULL, selected = NULL,
                                 options = list(placeholder = "Enter gene, eg: EGFR", 
                                                plugins = list('restore_on_backspace'))),
                  br(),
                  selectInput("plotTypeSel", 
                              label = h4("Plot type"), choices = ""),
                  br(),
                  actionButton(inputId = 'helpLink', label = 'Help', 
                               styleclass = "link", icon = "question-sign"),
#                   bsTooltip("dataset", "Select a dataset", "right", trigger="hover"),
#                   bsTooltip("gene", "Enter gene name", "right", trigger="hover"),
#                   bsTooltip("plotTypeSel", "Select one of the available plot", "right", trigger="hover"),
                  bsPopover(id = "helpLink", title =  "Help", 
                            content = "I am just trying to help here. I have to see if I can include a Rmd file ", 
                            trigger="hover", placement="bottom")
#                   uiOutput("help")
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
                                        checkboxInput("statTable", "Tukey's HSD", FALSE),
                                        checkboxInput("tTest", "Pairwise t tests", FALSE)))),
                             
                             plotOutput("plot"),
                             br(),
                             conditionalPanel(
                               condition = "input.statTable",
                               
                               strong("Tukey's Honest Significant Difference (HSD)"),
                               helpText("The table shows the difference between pairs, the 95% confidence interval and the p-value of the pairwise comparisons:"),
                               checkboxInput("stat", "Show the results in the plot", FALSE),
                               verbatimTextOutput("tukeyTest")),
                             conditionalPanel(
                               condition = "input.tTest",
                               br(),
                               strong("Pairwise t tests"),
                               helpText("Pairwise comparisons between group levels with corrections for multiple testing (p-values with Bonferroni correction):"),
                               verbatimTextOutput("pairwiseTtest")),
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
                                 max = 100),
                              
                               br(),
                               br(),
                               actionButton(inputId = 'downloadPlot', label = 'Download', 
                                            styleclass = "info", icon = "download-alt"))),
                    
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
                                                 conditionalPanel(
                                                   condition = "input.histologySurv == 'GBM'",
                                                   selectInput("subtypeSurv", strong("Subtype:"), 
                                                             c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP")))
                                          )
                                        ),
                                        conditionalPanel(
                                          condition = "input.histologySurv == 'GBM'",
                                          checkboxInput("gcimpSurv", "Exclude G-CIMP samples", FALSE)
                                        ),
                                        plotOutput("survPlot"),
                                        br(),
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
                                                 selectizeInput("gene1", h5("Gene 1"), "", options = list(plugins = list('restore_on_backspace')))
                                          ),
                                          column(4, 
                                                 selectizeInput("gene2", h5("Gene 2"), "", options = list(plugins = list('restore_on_backspace')))
                                          )),
                                        br(),
                                        fluidRow(
                                          column(3, 
                                                 selectInput("histologyCorr", h5("Histology:"), choices = "")
                                          ),
                                          column(3, 
                                                 conditionalPanel(
                                                   condition = "input.histologyCorr == 'GBM'",
                                                   selectInput("subtype", h5("Subtype (GBM):"), 
                                                             c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP")))
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
                             fluidRow(
                               splitLayout(
                                      uiOutput("piePlots"),
                                      uiOutput("survPlots")
                                      )),
#                              verbatimTextOutput("summary"),
                             br(),
                             downloadButton('downloadData', 'Download table'), 
                             br(),
                             br(),
                             dataTableOutput("table")),
                    
                    tabPanel("About", icon = icon("info-circle"),
                             includeMarkdown("tools/about.Rmd"))
                    
                  )
                )
  )
))
