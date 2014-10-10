# UI-elements for Tools tab
tabPanel(title = "Tools", icon = icon("gear"),
         tabsetPanel(
           
           tabPanel(title ="SubtypeME",
                    p(class = "lead","Classify tumor samples based on mRNA expression profiles"),
                    sidebarLayout(
                      sidebarPanel(width = 3,
                                   helpText("Upload a .csv file with samples in rows and genes expression in columns.",
                                            "The first column should contain the sample ID and should be named 'Sample'"),
                                   radioButtons(inputId = "tumorType", label = strong("Training set:"), 
                                                choices = c("GBM" = "gbm", "LGG" = "lgg")),
                                   fileInput(inputId = "upFile", label = "Choose CSV File",
                                             accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                   hr(),
                                   checkboxInput(inputId = "header", label = "Header", value = TRUE),
                                   radioButtons(inputId = "sep", label = "Separator",  
                                                choices = c(Comma = ",", Semicolon = ";", Tab="\t"), selected = ","),
                                   radioButtons(inputId = "quote", label = "Quote", 
                                                choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"'),
                                   hr(),
                                   br(),
                                   conditionalPanel(
                                     condition = "input.tabTools == 1",
                                     actionButton(inputId = "goSvm", label = "Submit SVM", styleclass = "primary"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 2",
                                     actionButton(inputId = "goGsva", label = "Submit ssGSEA", styleclass = "primary"),
                                     br()
                                   )
                      ),
                      mainPanel(
                        tabsetPanel(id = "tabTools",
                                    
                          tabPanel(title = "SVM", id = "SVM", value = 1,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="ajax-loader.gif"),
                                       hr(),
                                       p("Be patients, switching to another tab will crash GlioVis ...")
                                   ),
                                   p(class = "lead","Supported vector machine learning"),
                                   p(strong("File output format:"), "Download a .csv file containing the subtype call with the probability model."),  
                                   conditionalPanel(
                                     condition = "output.svm",
                                     downloadButton(outputId = "downloadSvm", label = "Download data")
                                   ),
                                   br(),
                                   dataTableOutput(outputId = "svm")
                          ),
                          
                          tabPanel(title = "ssGSEA", id = "GSVA", value = 2,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="ajax-loader.gif"),
                                       hr(),
                                       p("Be patients, switching to another tab will crash GlioVis ...")
                                   ),
                                   p(class = "lead","Single sample Gene Set Enrichment Analysis"),
                                   p(strong("File output format:"), "Download a .csv file containing the subtype call with the enrichment score."),
                                   conditionalPanel(
                                     condition = "output.gsva",
                                     downloadButton(outputId = "downloadGsva", label = "Download data")
                                   ),
                                   br(),
                                   dataTableOutput(outputId = "gsva")
                          )
                        )
                      )
                    )
           ),
           
           tabPanel(title = "CorrelateME", 
                    p(class = "lead","Correlate expression of a gene with all the genes in the dataset"),
                    br(),
                    sidebarLayout(
                      sidebarPanel(width = 3,
                                   selectInput(inputId = "datasetCor", label = h4("Dataset"),
                                               choices = datasets,
                                               selected = "TCGA GBM", selectize = TRUE),
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
                                   # it triggers after the gene is selected
                                   # probably need to use isolate()
                                   br(),
                                   conditionalPanel(
                                     condition = "output.corrData",
                                     br(),
                                     downloadButton(outputId = "downloadCorrData", label = "Download data")
                                   )
                      ),
                      mainPanel(
                        column(width = 8,
                               dataTableOutput(outputId = "corrData")
                        )
                      )
                    )
           )
         )
)