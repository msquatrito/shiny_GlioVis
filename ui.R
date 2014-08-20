library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)
library (shinysky)
library (shinyBS)
library (dplyr)
# library (kernlab)

shinyUI(
  
  fluidPage(
    
    # progressInit() must be called somewhere in the UI in order for the progress UI to actually appear
    progressInit(),
    
    # link to the style.css file. Currently need only for the positioning of the progress bar
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
    
    sidebarLayout(
      
      sidebarPanel(
        img(src = "GlioVis_logo.jpg", height = 90, width = 270),
        br(),
        br(),
        selectInput(inputId = "dataset", label = h4("Dataset"),
                    choices = c("TCGA GBM", "TCGA Lgg", "Rembrandt","Gravendeel", "Phillips", "Murat", "Freije"),
                    selected = "TCGA GBM", selectize = FALSE),
        br(),
        selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL, selected = NULL,
                       options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),
        br(),
        selectInput(inputId = "plotTypeSel", label = h4("Plot type"), choices = "", selectize = FALSE),
        br(),
        #         actionButton(inputId = 'helpLink', label = 'Help', styleclass = "link", icon = "question-sign"),
        bsTooltip("dataset", "Choose a dataset", "right", trigger="hover"),
        bsTooltip("gene", "Enter gene name", "right", trigger="hover"),
        bsTooltip("plotTypeSel", "Select one of the available plot for the specified dataset", "right", trigger="hover"),
        bsPopover(id = "helpLink", title =  "Help", 
                  content = "I am just trying to help here. I have to see if I can include a Rmd file ", 
                  trigger="hover", placement="bottom")
        #         helpPopup(title = "Help me pleaseeeeee", content = includeMarkdown("tools/help.Rmd"), 
        #                   placement = "bottom", trigger = "click")
        #         uiOutput("help")
      ),
      
      
      mainPanel(
        
        h3(textOutput(outputId = "caption")),
        
        tabsetPanel(
          
          tabPanel(title = "Home", icon = icon("home"),
                   #                    includeMarkdown("tools/home.Rmd"),
                   #                    hr(),
                   helpText("HOME TAB UNDER CONSTRUCTION"),
                   helpText("In this tab I would add a wellcome message, some info regarding available datasets and  plots"),
                   hr(),
                   h4("Datasets"),                   
                   img(src = "dataset.jpg", height = 200, width = 800),
                   hr(),
                   h4("Funding"),
                   img(src = "seve.jpg", height = 135, width = 405)
                   
          ),
          
          tabPanel(title = "Plot", icon = icon("bar-chart-o"), id = "plots",
                   sidebarPanel(
                     h5("Plot options:"),
                     checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                     checkboxInput(inputId = "colStrip", label = "Color stripchart", value = FALSE),
                     checkboxInput(inputId = "colBox", label = "Color box", value = FALSE),
                     checkboxInput(inputId = "bw", label = "White background", value = FALSE),
                     hr(),
                     h5("Statistic:"),
                     checkboxInput(inputId = "statTable", label = "Tukey's HSD", value = FALSE),
                     checkboxInput(inputId = "tTest", label = "Pairwise t tests", value = FALSE),
                     hr(),
                     selectInput(inputId = "downloadPlotFileType", label = strong("Select download file type"),
                                 choices = list("PDF"  = "pdf", "BMP"  = "bmp", "JPEG" = "jpeg", "PNG"  = "png")
                     ),
                     # Allow the user to set the height and width of the plot download.
                     helpText(strong("Set download image dimensions"),
                              "(units are inches for PDF, pixels for all other formats)"),
                     numericInput(inputId = "downloadPlotHeight", label = "Height (inches)", value = 7, min = 1, max = 100),
                     numericInput(inputId = "downloadPlotWidth", label = "Width (inches)", value = 7, min = 1, max = 100),
                     br(),
                     br(),
                     downloadButton(outputId = "downloadPlot", label = "Download")
                   ),
                   mainPanel(
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
                     )
                   )
          ),
          
          tabPanel(title = "Survival", icon = icon("user-md"),
                   tabsetPanel(
                     tabPanel(title = "Km plot",
                              sidebarPanel( 
                                selectInput(inputId = "histologySurv", label = h5("Histology:"), choices = ""),
                                conditionalPanel(
                                  condition = "input.histologySurv == 'GBM'",
                                  checkboxInput(inputId = "gcimpSurv", label = "Exclude G-CIMP samples", value = FALSE)
                                ),
                                conditionalPanel(
                                  condition = "input.histologySurv == 'GBM'",
                                  selectInput(inputId = "subtypeSurv", label = h5("Subtype:"), 
                                              choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                                ),
                                hr(),
                                selectInput(inputId = "cutoff", label = h5("Select cutoff:"), 
                                            choices = c("median", "lower quartile", "upper quartile", "quartiles")),
                                hr(),
                                downloadButton(outputId = "downloadsurvPlot", label = "Download")
                              ),
                              mainPanel(
                                plotOutput(outputId = "survPlot")
                              )
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
                              sidebarPanel(
                                selectizeInput(inputId = "gene1", label = h5("Gene 1"), choices ="", 
                                               options = list(plugins = list('restore_on_backspace'))),
                                selectizeInput(inputId = "gene2", label = h5("Gene 2"), choices ="", 
                                               options = list(plugins = list('restore_on_backspace'))),
                                hr(),
                                selectInput(inputId = "histologyCorr", label = h5("Histology:"), choices = ""),
                                conditionalPanel(
                                  condition = "input.histologyCorr == 'GBM'",
                                  selectInput(inputId = "subtype", label = h5("Subtype (GBM):"), 
                                              choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                                ),
                                hr(),
                                h5("Statistic:"),
                                checkboxInput(inputId = "statCorr", label = "Pearson's correlation", value = FALSE),
                                hr(),
                                radioButtons(inputId = "colorBy", label = h5("Color by:"), 
                                             choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                                hr(),
                                radioButtons(inputId = "separateBy", label = h5("Separate by:"), 
                                             choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                                hr(),
                                downloadButton(outputId = "downloadcorrPlot", label = "Download plot")
                              ),
                              
                              mainPanel(
                                plotOutput(outputId = "corrPlot"), 
                                br(),
                                conditionalPanel(
                                  condition = "input.statCorr",
                                  verbatimTextOutput(outputId = "corrTest")
                                )
                              )
                     ),
                     
                     tabPanel(title = "Correlation data", 
                              sidebarPanel( 
                                selectInput(inputId = "histologyCorrTable", label = "Histology:", choices = ""),
                                hr(),
                                selectInput(inputId = "cor", label = "Correlation:", choices = c("All", "Positive", "Negative")),
                                hr(),
                                selectInput(inputId = "sign", label = "Signficance:", choices = c(0.05, 0.01)),
                                hr(),
                                downloadButton(outputId = "downloadCorrData", label = "Download data")
                              ),
                              mainPanel(
                                dataTableOutput(outputId = "corrData")
                              )
                     )
                   )
          ),
          
          tabPanel(title = "Data", icon = icon("table"),
                   tabsetPanel(
                     tabPanel(title = "Table",
                              downloadButton(outputId = "downloadData", label = "Download table"), 
                              br(),
                              br(),
                              dataTableOutput(outputId = "table")
                     ),
                     tabPanel(title = "Summary plots",
                              splitLayout(
                                uiOutput(outputId = "piePlots", inline = TRUE), # inline = TRUE not working
                                uiOutput(outputId = "survPlots", inline = TRUE)
                              )
                     )
                   )
          ),
          
          tabPanel(title = "Tools", icon = icon("gear"),
                   tabsetPanel(
                     tabPanel(title = "SubtypeME",
                              helpText(HTML("Classify tumor samples based on mRNA expression profiles using Supported Vector Machine Learning")),
                              helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                              helpText(HTML("<b>File input format: </b> Upload a .csv file with samples in rows and genes expression in columns.<br>
                                            The first column should contain the sample ID and should be named 'Sample'.<br>")),
                              sidebarPanel(
                                radioButtons(inputId = "svm", label = NULL, 
                                             choices = c("GBM" = "gbm", "LGG (NOT ACTIVE YET)" = "lgg")),
                                fileInput(inputId = "upFile", label = "Choose CSV File",
                                          accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                hr(),
                                checkboxInput(inputId = "header", label = "Header", value = TRUE),
                                radioButtons(inputId = "sep", label = "Separator",  
                                             choices = c(Comma = ",", Semicolon = ";", Tab="\t"), selected = ","),
                                radioButtons(inputId = "quote", label = "Quote", 
                                             choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"')
                              ),
                              mainPanel(
                                tableOutput(outputId = "svm")
                              )
                     )
                   )
          ),
          
          tabPanel(title = "About", icon = icon("info-circle"), 
                   includeMarkdown("tools/about.Rmd"))
          
        )
        )
    )
  )
  )
