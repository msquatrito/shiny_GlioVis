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
                                     actionButton(inputId = "goKnn", label = "Submit K-NN", styleclass = "primary"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 3",
                                     actionButton(inputId = "goGsva", label = "Submit ssGSEA", styleclass = "primary"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 4",
                                     actionButton(inputId = "goSub3", label = "Submit 3-Way", styleclass = "primary"),
                                     br()
                                   )
                      ),
                      mainPanel(
                        tabsetPanel(id = "tabTools",
                                    
                          tabPanel(title = "SVM", id = "SVM", value = 1,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="Rotating_brain.gif"),
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
                          
                          tabPanel(title = "K-NN", id = "K-NN", value = 2,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="Rotating_brain.gif"),
                                       hr(),
                                       p("Be patients, switching to another tab will crash GlioVis ...")
                                   ),
                                   p(class = "lead","K-nearest neighbors prediction"),
                                   p(strong("File output format:"), "Download a .csv file containing the subtype call with the probability score."),
                                   conditionalPanel(
                                     condition = "output.knn",
                                     downloadButton(outputId = "downloadKnn", label = "Download data")
                                   ),
                                   br(),
                                   dataTableOutput(outputId = "knn")
                          ),
                          
                          tabPanel(title = "ssGSEA", id = "GSVA", value = 3,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="Rotating_brain.gif"),
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
                          ),
                          
                          tabPanel(title = "3-Way", id = "3-Way", value = 4,
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="Rotating_brain.gif"),
                                       hr(),
                                       p("Be patients, it could take a couple of minutes ...")
                                   ),
                                   p(class = "lead","Generate and compare subtype calls by SVM, K-NN and ssGSEA"),
                                   conditionalPanel(
                                     condition = "output.sub3",
                                     p(strong("svm_vs_knn"),textOutput(outputId = "call.identity",inline = TRUE)),
                                     tableOutput(outputId = "sub3Summary.1"),
                                     br(),
                                     fluidRow(
                                       column(width = 5,
                                              p(strong("svm_vs_gsea")),
                                              tableOutput(outputId = "sub3Summary.2")
                                       ),
                                       column(width = 5, offset = 1,
                                              p(strong("knn_vs_gsea")),
                                              tableOutput(outputId = "sub3Summary.3")
                                       )
                                     )
                                   ),
                                   br(),
                                   conditionalPanel(
                                     condition = "output.sub3",
                                     downloadButton(outputId = "downloadSub3", label = "Download data")
                                   ),
                                   br(),
                                   dataTableOutput(outputId = "sub3")
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
                                   div(class = "busy",  
                                       p("Calculating, please wait"),
                                       img(src="Rotating_brain.gif") 
                                   ),
                                   hr(),
                                   selectInput(inputId = "cor", label = h5("Correlation:"), choices = c("All", "Positive", "Negative")),
                                   br(),
                                   selectInput(inputId = "sign", label = h5("Signficance:"), choices = c(0.05, 0.01)),
                                   br(),
                                   selectInput(inputId = "corrMethod",label = h5("Method:"), choices = c("Pearson", "Spearman")),
                                   hr(),
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