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


shinyUI(
  
  navbarPage(title = "GlioVis", windowTitle = "GlioVis - Visualization Tools for Glioma Datasets", fluid = TRUE,
             
             tabPanel(title = "Home", icon = icon("home"),
                      div(class = "outer", 
#                       tags$head(
#                         tags$style(type="text/css",'body {
#                         background-image: url("gradient_inv.jpg");
#                         background-repeat: repeat-x;
#                       }')
#                       ),
              div(class = "span9",               
                      img(src = "GlioVis_logo_white.jpg", width = 750),
                      hr(),
                      p(class="lead","Welcome to", strong("GlioVis"),": a user friendly web based data visualization and analysis application for exploring Glioma datasets."),
                      br(),
                      strong("Available datasets:"),                   
                      source("tools/datsets.table.R",local = TRUE)$value,
                      br(),br(),
                      p(strong("Please cite:"), a("Bowman R. and Squatrito M.", href="", target="_blank"), " (manuscript in preparation) when publishing results based on GlioVis."),
                      hr(),
                      tags$blockquote(# class="pull-right",
                        tags$p("No great discovery was ever made without a bold guess."), 
                        tags$small("Isaac Newton"))
              )
                      )
                      
             ),
             
             tabPanel(title = "Explore", icon = icon("picture-o"), id = "explore",
                      
                      # progressInit() must be called somewhere in the UI in order for the progress UI to actually appear
                      progressInit(),
                      
                      # link to the style.css file. Currently need only for the positioning of the progress bar
                      tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                      
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     img(src = "GlioVis_logo_gray.jpg"),
                                     br(),
                                     br(),
                                     selectInput(inputId = "dataset", label = h4("Dataset"),
                                                 choices = c("TCGA GBM", "TCGA Lgg", "Rembrandt","Gravendeel", "Phillips", "Murat", "Freije"),
                                                 selected = "TCGA GBM", selectize = FALSE),
                                     br(),
                                     selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL, selected = NULL,
                                                    options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),

                                     conditionalPanel(
                                       condition = "input.tab1 == 1",
                                       br(),
                                       selectInput(inputId = "plotTypeSel", label = h4("Plot type"), choices = "", selectize = FALSE),
                                       hr(),
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
                                     conditionalPanel(
                                       condition = "input.tab1 == 2",
                                       br(),
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
                                     conditionalPanel(
                                       condition = "input.tab1 == 3",
                                       br(),
                                       selectizeInput(inputId = "gene2", label = h4("Gene 2"), choices ="", 
                                                      options = list(placeholder = "Enter gene 2, eg: SOCS2", plugins = list('restore_on_backspace'))),
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
                          
                          tabsetPanel(id = "tab1",
                                      
                            tabPanel(title = "Box plots", icon = icon("bar-chart-o"), id = "plots", value = 1,
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
                            ),
                            
                            tabPanel(title = "Survival", icon = icon("user-md"), value = 2,
                                     tabsetPanel(id = "tabSurv",
                                       tabPanel(title = "Km plot",  value = "km",
                                                plotOutput(outputId = "survPlot", width = 500 , height = 400)
                                       ),
                                       tabPanel(title = "HR plot", value = "hr",
                                                wellPanel(
                                                  helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                                                  helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                              to update the survival plot. The blue line represents the current selection.")),
                                                  br(),
                                                  checkboxInput(inputId = "quantile", label = "Show quantiles", value = TRUE),
                                                  checkboxInput(inputId = "gcimp", label = "Exclude G-CIMP samples", value = FALSE),
                                                  plotOutput(outputId = "hazardPlot", clickId = "densityClick", width = 500 , height = 400)
                                                ),
                                                plotOutput(outputId = "kmPlot", width = 500 , height = 400),
                                                downloadButton(outputId = "downloadkmPlot", label = "Download")
                                                
                                       )
                                     )
                            ),
                                                       
                            tabPanel(title = "Correlations", icon = icon("list-alt"),  value = 3,
                                     plotOutput(outputId = "corrPlot"), 
                                     br(),
                                     conditionalPanel(
                                       condition = "input.statCorr",
                                       verbatimTextOutput(outputId = "corrTest")
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
                                                  uiOutput(outputId = "piePlots", inline = TRUE),# inline = TRUE not working
                                                  uiOutput(outputId = "survPlots", inline = TRUE)
                                                )
                                       )
                                     )
                            )                  
                          )
                        )
                      )
             ),
             
             tabPanel(title = "Tools", icon = icon("gear"),
                      tabsetPanel(
                        tabPanel(title = "SubtypeME",
                                 p(class = "lead","Classify tumor samples based on mRNA expression profiles using Supported Vector Machine Learning"),
                                 p(strong("IMPORTANT:"),"Currently active only for GBM samples. Check the 'About' tab for more information."),
                                 p(strong("File input format:"), "Upload a .csv file with samples in rows and genes expression in columns.The first column should contain the sample ID and should be named 'Sample'"),
                                 p(strong("File output format:"), "Download a .csv file containing the subtype call with the probability model."),
                                 br(),
                                 sidebarLayout(
                                   sidebarPanel(width = 3,
                                                radioButtons(inputId = "svm", label = NULL, 
                                                             choices = c("GBM" = "gbm", "LGG (NOT ACTIVE YET)" = "lgg")),
                                                fileInput(inputId = "upFile", label = "Choose CSV File",
                                                          accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                                hr(),
                                                checkboxInput(inputId = "header", label = "Header", value = TRUE),
                                                radioButtons(inputId = "sep", label = "Separator",  
                                                             choices = c(Comma = ",", Semicolon = ";", Tab="\t"), selected = ","),
                                                radioButtons(inputId = "quote", label = "Quote", 
                                                             choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"'),
                                                hr(),
                                                conditionalPanel(
                                                  condition = "output.svm",
                                                  downloadButton(outputId = "downloadSvm", label = "Download data"),
                                                  br()
                                                )
                                   ),
                                   mainPanel(
                                     dataTableOutput(outputId = "svm")
                                   )
                                 )
                        ),
                        tabPanel(title = "CorrelateME", 
                                 p(class = "lead","Correlate expression of a gene with all the genes in the dataset"),
                                 helpText("Calculating gene correlations will take a couple of minutes. Please DON'T switch tabs while computing"),
                                 br(),
                                 sidebarLayout(
                                   sidebarPanel(width = 3,
                                                selectInput(inputId = "datasetCor", label = h4("Dataset"),
                                                            choices = c("TCGA GBM", "TCGA Lgg", "Rembrandt","Gravendeel", "Phillips", "Murat", "Freije"),
                                                            selected = "TCGA GBM", selectize = FALSE),
                                                br(),
                                                selectInput(inputId = "histologyCorrTable", label = h4("Histology:"), choices = ""),
                                                br(),
                                                selectizeInput(inputId = "geneCor", label = h4("Gene"), choices = NULL, selected = NULL,
                                                               options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),
                                                hr(),
                                                selectInput(inputId = "cor", label = h5("Correlation:"), choices = c("All", "Positive", "Negative")),
                                                br(),
                                                selectInput(inputId = "sign", label = h5("Signficance:"), choices = c(0.05, 0.01)),
                                                br(),
                                                actionButton(inputId = "goCor", label = "Get data", styleclass = "primary"), # Not working correctly, the second time 
                                                                                                                            # it triggers after the gene si selected
                                                                                                                            # probably need to use isolate()
                                                br(),
                                                conditionalPanel(
                                                  condition = "output.corrData",
                                                  downloadButton(outputId = "downloadCorrData", label = "Download data")
                                                )
                                   ),
                                   mainPanel(
                                     dataTableOutput(outputId = "corrData")
                                   )
                                 )
                        )
                      )
             ),
             
             tabPanel(title = "About", icon = icon("info-circle"), 
                      includeMarkdown("tools/about.Rmd")
                      )
  )
)
