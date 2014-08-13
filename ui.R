library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)
library (shinysky)
library (shinyBS)
library (kernlab)
library (dplyr)


shinyUI(
  
  fluidPage(
    
    # progressInit() must be called somewhere in the UI in order for the progress UI to actually appear
    progressInit(),

    sidebarLayout(
      
      sidebarPanel(
        
        img(src = "GlioVis_logo.jpg", height = 90, width = 270),
        br(),
        br(),
        selectInput(inputId = "dataset", label = h4("Dataset"),
                    choices = c("TCGA GBM", "TCGA Lgg", "Rembrandt","Gravendeel", "Phillips", "Murat", "Freije"),
                    selected = "TCGA GBM"),
        br(),
        selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL, selected = NULL,
                       options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),
        br(),
        selectInput(inputId = "plotTypeSel", label = h4("Plot type"), choices = ""),
        br(),
#         actionButton(inputId = 'helpLink', label = 'Help', styleclass = "link", icon = "question-sign"),
#                   bsTooltip("dataset", "Select a dataset", "right", trigger="hover"),
#                   bsTooltip("gene", "Enter gene name", "right", trigger="hover"),
#                   bsTooltip("plotTypeSel", "Select one of the available plot", "right", trigger="hover"),
#         bsPopover(id = "helpLink", title =  "Help", 
#                   content = "I am just trying to help here. I have to see if I can include a Rmd file ", 
#                   trigger="hover", placement="bottom")
#                   helpPopup(title = "Help me pleaseeeeee", content = includeMarkdown("tools/help.Rmd"), 
#                             placement = "bottom", trigger = "click")
        uiOutput("help")
      ),
      
      
      mainPanel(
        
        h3(textOutput(outputId = "caption")),
        
        tabsetPanel(
          
          tabPanel(title = "Plot", icon = icon("bar-chart-o"), id = "plots",
                   fluidRow(
                     column(width = 4, 
                            wellPanel(
                              strong("Plot options:"),
                              checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                              checkboxInput(inputId = "colStrip", label = "Color stripchart", value = FALSE),
                              checkboxInput(inputId = "colBox", label = "Color box", value = FALSE)
                            )
                     ),    
                     column(width = 4, 
                            wellPanel(
                              strong("Statistic:"),
                              checkboxInput(inputId = "statTable", label = "Tukey's HSD", value = FALSE),
                              checkboxInput(inputId = "tTest", label = "Pairwise t tests", value = FALSE)
                            )
                     )
                   ),
                   plotOutput(outputId = "plot"),
                   br(),
                   conditionalPanel(
                     condition = "input.statTable",
                     strong("Tukey's Honest Significant Difference (HSD)"),
                     helpText("The table shows the difference between pairs, the 95% confidence interval and the p-value of the pairwise comparisons:"),
                     checkboxInput(inputId = "stat", label = "Show the results in the plot", value = FALSE),
                     verbatimTextOutput(outputId = "tukeyTest")
                   ),
                   conditionalPanel(
                     condition = "input.tTest",
                     br(),
                     strong("Pairwise t tests"),
                     helpText("Pairwise comparisons between group levels with corrections for multiple testing (p-values with Bonferroni correction):"),
                     verbatimTextOutput(outputId = "pairwiseTtest")
                   ),
                   br(),
                   wellPanel(
                     selectInput(inputId = "downloadPlotFileType", label   = h5("Select download file type"),
                                 choices = list("PDF"  = "pdf", "BMP"  = "bmp", "JPEG" = "jpeg", "PNG"  = "png")
                     ),
                     # Allow the user to set the height and width of the plot download.
                     h5(HTML("Set download image dimensions<br>(units are inches for PDF, pixels for all other formats)")),
                     numericInput(inputId = "downloadPlotHeight", label = "Height (inches)", value = 7, min = 1, max = 100),
                     numericInput(inputId = "downloadPlotWidth", label = "Width (inches)", value = 7, min = 1, max = 100),
                     br(),
                     br(),
                     downloadButton(outputId = "downloadPlot", label = "Download")
                   )
          ),
          
          tabPanel(title = "Survival", icon = icon("user-md"),
                   tabsetPanel(
                     tabPanel(title = "Km plot",
                              fluidRow(
                                column(width = 4, 
                                       selectInput(inputId = "histologySurv", label = strong("Histology:"), choices = "")
                                ),
                                column(width = 4,
                                       selectInput(inputId = "cutoff", label = strong("Select cutoff:"), 
                                                   choices = c("median", "lower quartile", "upper quartile", "quartiles"))
                                ),
                                column(width = 4, 
                                       conditionalPanel(
                                         condition = "input.histologySurv == 'GBM'",
                                         selectInput(inputId = "subtypeSurv", label = strong("Subtype:"), 
                                                     choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                                       )
                                )
                              ),
                              conditionalPanel(
                                condition = "input.histologySurv == 'GBM'",
                                checkboxInput(inputId = "gcimpSurv", label = "Exclude G-CIMP samples", value = FALSE)
                              ),
                              plotOutput(outputId = "survPlot"),
                              br(),
                              downloadButton(outputId = "downloadsurvPlot", label = "Download")
                     ),
                     
                     tabPanel(title = "HR plot", 
                              wellPanel(
                                helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                                helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                              to update the survival plot. The blue line represents the current selection.")),
                                br(),
                                checkboxInput(inputId = "quantile", label = "Show quantiles", value = TRUE),
                                checkboxInput(inputId = "gcimp", label = "Exclude G-CIMP samples", value = FALSE),
                                plotOutput(outputId = "hazardPlot", clickId = "densityClick")
                              ),
                              plotOutput(outputId = "kmPlot"),
                              downloadButton(outputId = "downloadkmPlot", label = "Download")
                     )
                   )
          ),
          
          
          tabPanel(title = "Correlations", icon = icon("list-alt"),  
                   tabsetPanel(
                     tabPanel(title = "Correlation plot",  
                              fluidRow(
                                column(width = 4,
                                       selectizeInput("gene1", h5("Gene 1"), "", options = list(plugins = list('restore_on_backspace')))
                                ),
                                column(width = 4, 
                                       selectizeInput("gene2", h5("Gene 2"), "", options = list(plugins = list('restore_on_backspace')))
                                )
                              ),
                              br(),
                              fluidRow(
                                column(width = 3, 
                                       selectInput("histologyCorr", h5("Histology:"), choices = "")
                                ),
                                column(width = 3, 
                                       conditionalPanel(
                                         condition = "input.histologyCorr == 'GBM'",
                                         selectInput(inputId = "subtype", label = h5("Subtype (GBM):"), 
                                                     choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP")))
                                ),
                                column(width = 3, 
                                       inputPanel(
                                         radioButtons(inputId = "colorBy", label = strong("Color by:"), 
                                                      choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"))
                                       )
                                ),
                                column(width = 3, 
                                       inputPanel(
                                         radioButtons(inputId = "separateBy", label = strong("Separate by:"), 
                                                      choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"))
                                       )
                                )
                              ),
                              br(),
                              plotOutput(outputId = "corrPlot"),
                              downloadButton(outputId = "downloadcorrPlot", label = "Download plot"), 
                              br(),
                              br(),
                              verbatimTextOutput(outputId = "corrTest")
                     ),
                     
                     tabPanel(title = "Correlation data", 
                              fluidRow(
                                column(width = 4, 
                                       selectInput(inputId = "histologyCorrTable", label = "Histology:", choices = "")
                                ),
                                column(width = 4, 
                                       selectInput(inputId = "cor", label = "Correlation:", choices = c("All", "Positive", "Negative"))
                                ),
                                column(width = 4, 
                                       selectInput(inputId = "sign", label = "Signficance:", choices = c(0.05, 0.01))
                                )
                              ),
                              br(),
                              downloadButton(outputId = "downloadCorrData", label = "Download data"),
                              dataTableOutput(outputId = "corrData")
                     )
                   )
          ),
          
          tabPanel(title = "Data", icon = icon("table"),
                   fluidRow(
                     splitLayout(
                       uiOutput(outputId = "piePlots", inline = TRUE), # inline = TRUE not working
                       uiOutput(outputId = "survPlots", inline = TRUE)
                     )
                   ),
                   #                              verbatimTextOutput("summary"),
                   br(),
                   downloadButton(outputId = "downloadData", label = "Download table"), 
                   br(),
                   br(),
                   dataTableOutput(outputId = "table")
          ),
          
          tabPanel(title = "Tools", icon = icon("gear"),
                   tabsetPanel(
                     tabPanel(title = "SubtypeME",
                              helpText(HTML("Classify GBM samples based on mRNA expression profiles using Supported Vector Machine Learning")),
#                               helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                              helpText(HTML("<b>File input format: </b> Upload a .csv file with samples in rows and genes expression in columns.<br>
                                            The first column should contain the sample ID and should be named 'Sample'.<br>")),
                              fluidRow(
                                column(width = 4,
                                       wellPanel(width = 400,

                                        fileInput('upFile', 'Choose CSV File',
                                                  accept=c('text/csv', 
                                                           'text/comma-separated-values,text/plain', 
                                                           '.csv')),
                                        tags$hr(),
                                        checkboxInput('header', 'Header', TRUE),
                                        radioButtons('sep', 'Separator',
                                                     c(Comma=',',
                                                       Semicolon=';',
                                                       Tab='\t'),
                                                     ','),
                                        radioButtons('quote', 'Quote',
                                                     c(None='',
                                                       'Double Quote'='"',
                                                       'Single Quote'="'"),
                                                     '"'))
                                       ),
                                column(width = 8,
                                       tableOutput('svm'))
                              )
                     )
                   )
          ),
          
          tabPanel(title = "About", icon = icon("info-circle"), includeMarkdown("tools/about.Rmd"))
          
        )
      )
    )
  )
)
