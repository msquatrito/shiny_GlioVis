# UI-elements for Tools tab
tabPanel(title = "Tools", icon = icon("gear"),
         tabsetPanel(
           
           tabPanel(title = "CorrelateME", 
                    p(class = "lead","Correlate expression of a gene with all the genes in the dataset"),
                    sidebarLayout(
                      sidebarPanel(width = 3,
                                   selectInput(inputId = "datasetCor", label = h4("Dataset"),
                                               choices = datasets,
                                               selected = "TCGA GBM", selectize = TRUE),
                                   selectInput(inputId = "histologyCorrTable", label = h4("Histology:"), choices = ""),
                                   selectizeInput(inputId = "geneCor", label = h4("Gene"), choices = NULL, selected = NULL,
                                                  options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),
                                   hr(),
                                   selectInput(inputId = "cor", label = h5("Correlation:"), choices = c("All", "Positive", "Negative")),
                                   conditionalPanel(
                                     condition = "input.cor == 'All'",
                                     sliderInput("range", "Range:",min = -1, max = 1, value = c(-0.3,0.3),step = 0.1,round = FALSE)
                                   ),
                                   radioButtons(inputId = "sign", label = h5("Signficance:"), choices = c(0.05, 0.01),inline = TRUE),
                                   radioButtons(inputId = "corrMethod",label = h5("Method:"), choices = c("Pearson", "Spearman"),inline = TRUE),
                                   hr(),
                                   conditionalPanel(
                                     condition = "output.corrData",
                                     br(),
                                     downloadButton(outputId = "downloadCorrData", label = "Download data", class= "btn-primary")
                                   )
                      ),
                      mainPanel(
                        div(class = "busy",  
                            p("Calculating, please wait"),
                            img(src="Rotating_brain.gif") 
                        ),
                        splitLayout(cellWidths = c("60%", "40%"),
                                    dataTableOutput(outputId = "corrData"),
                                    plotOutput(outputId = "corrDataPlot")
                        )
                        #                         column(width = 8,
                        #                                dataTableOutput(outputId = "corrData")
                        #                         ),
                        #                         column(width = 4,
                        #                                br(),
                        #                                br(),
                        #                                plotOutput(outputId = "corrDataPlot")
                        #                         )
                      )
                    )
           ),
           
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
                                     condition = "input.tabTools == 1 & output.finishedUploading",
                                     actionButton(inputId = "goSvm", label = "Submit SVM", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 2 & output.finishedUploading",
                                     actionButton(inputId = "goKnn", label = "Submit K-NN", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 3 & output.finishedUploading",
                                     actionButton(inputId = "goGsva", label = "Submit ssGSEA", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabTools == 4 & output.finishedUploading",
                                     actionButton(inputId = "goSub3", label = "Submit 3-Way", class= "btn-success"),
                                     br()
                                   )
                      ),
                      mainPanel(
                        tabsetPanel(id = "tabTools",
                                    
                                    tabPanel(title = "SVM", id = "SVM", value = 1,
                                             p(class = "lead","Supported vector machine learning"),
                                             p(strong("File output format:"), "Subtype calls with the probability model."),  
                                             conditionalPanel(
                                               condition = "output.svm",
                                               downloadButton(outputId = "downloadSvm", label = "Download data", class= "btn-primary")
                                             ),
                                             br(),
                                             div(class = "busy",  
                                                 p("Calculating, please wait"),
                                                 img(src="Rotating_brain.gif"),
                                                 hr(),
                                                 p("Be patients, switching to another tab will crash GlioVis ...")
                                             ),
                                             dataTableOutput(outputId = "svm")
                                    ),
                                    
                                    tabPanel(title = "K-NN", id = "K-NN", value = 2,
                                             
                                             p(class = "lead","K-nearest neighbors prediction"),
                                             p(strong("File output format:"), "Subtype calls with the probability score."),
                                             conditionalPanel(
                                               condition = "output.knn",
                                               downloadButton(outputId = "downloadKnn", label = "Download data", class= "btn-primary")
                                             ),
                                             br(),
                                             div(class = "busy",  
                                                 p("Calculating, please wait"),
                                                 img(src="Rotating_brain.gif"),
                                                 hr(),
                                                 p("Be patients, switching to another tab will crash GlioVis ...")
                                             ),
                                             dataTableOutput(outputId = "knn")
                                    ),
                                    
                                    tabPanel(title = "ssGSEA", id = "GSVA", value = 3,
                                             p(class = "lead","Single sample Gene Set Enrichment Analysis"),
                                             p(strong("File output format:"), "Subtype calls with the enrichment score."),
                                             conditionalPanel(
                                               condition = "output.gsva",
                                               downloadButton(outputId = "downloadGsva", label = "Download data", class= "btn-primary")
                                             ),
                                             br(),
                                             div(class = "busy",  
                                                 p("Calculating, please wait"),
                                                 img(src="Rotating_brain.gif"),
                                                 hr(),
                                                 p("Be patients, switching to another tab will crash GlioVis ...")
                                             ),
                                             dataTableOutput(outputId = "gsva")
                                    ),
                                    
                                    tabPanel(title = "3-Way", id = "3-Way", value = 4,
                                             p(class = "lead","Generate and compare subtype calls by SVM, K-NN and ssGSEA"),
                                             br(),
                                             conditionalPanel(
                                               condition = "output.sub3",
                                               downloadButton(outputId = "downloadSub3", label = "Download data", class= "btn-primary")
                                             ),
                                             br(),
                                             div(class = "busy",  
                                                 p("Calculating, please wait"),
                                                 img(src="Rotating_brain.gif"),
                                                 hr(),
                                                 p("Be patients, switching to another tab will crash GlioVis ...")
                                             ),
                                             column(width = 9,
                                                    dataTableOutput(outputId = "sub3")
                                             )
                                    )
                        )
                      )
                    )
           )
         )
)