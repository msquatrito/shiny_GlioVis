# UI-elements for Tools tab
tabPanel(title = "Tools", icon = icon("gear"),
         tabsetPanel(
           tabPanel(title = "CorrelateME", 
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
                                   radioButtons(inputId = "corrMethod",label = h5("Method:"), choices = c("Pearson", "Spearman"),inline = TRUE)
                      ),
                      mainPanel(
                        p(class = "lead","Correlate expression of a gene with all the genes in the dataset"),
                        div(class = "busy",  
                            p("Calculating, please wait"),
                            img(src="Rotating_brain.gif") 
                        ),
                        splitLayout(cellWidths = c("60%", "40%"),
                                    DT::dataTableOutput(outputId = "corrData"),
                                    plotOutput(outputId = "corrDataPlot")
                        )
                      )
                    )
           ),
           
           tabPanel(title ="SubtypeME",
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
                                     condition = "input.tabClassify == 'SVM' & output.finishedUploading",
                                     actionButton(inputId = "goSvm", label = "Submit SVM", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabClassify == 'K-NN' & output.finishedUploading",
                                     actionButton(inputId = "goKnn", label = "Submit K-NN", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabClassify == 'ssGSEA' & output.finishedUploading",
                                     actionButton(inputId = "goGsva", label = "Submit ssGSEA", class= "btn-success"),
                                     br()
                                   ),
                                   conditionalPanel(
                                     condition = "input.tabClassify == '3-Way' & output.finishedUploading",
                                     actionButton(inputId = "goSub3", label = "Submit 3-Way", class= "btn-success"),
                                     br()
                                   )
                      ),
                      mainPanel(
                        div(style = "width: 100%; overflow: hidden;",
                            div(style = "width: 550px; float: left;",p(class = "lead","Classify tumor samples based on mRNA expression profiles")),
                            div(style = "margin-left: 550px;", helpModal(modal_title ="SubtypeMe", link = "helpClassify", help_file = includeMarkdown("tools/help_classify.Rmd")))
                        ),
                        
                        tabsetPanel(id = "tabClassify",   
                                    
                                    tabPanel(title = "SVM", id = "SVM", 
                                             p(class = "lead","Supported vector machine learning"),
                                             p(strong("Output format:"), "Subtype calls with the probability model."),  
                                             br(),
                                             column(width = 9,
                                                    div(class = "busy",  
                                                        p("Calculating, please wait"),
                                                        img(src="Rotating_brain.gif"),
                                                        hr(),
                                                        p("Be patients, switching to another tab will crash GlioVis ...")
                                                    ),
                                                    DT::dataTableOutput(outputId = "svm")
                                             )
                                    ),
                                    
                                    tabPanel(title = "K-NN", id = "K-NN",                                            
                                             p(class = "lead","K-nearest neighbors prediction"),
                                             p(strong("Output format:"), "Subtype calls with the probability score."),
                                             br(),
                                             column(width = 8,
                                                    div(class = "busy",  
                                                        p("Calculating, please wait"),
                                                        img(src="Rotating_brain.gif"),
                                                        hr(),
                                                        p("Be patients, switching to another tab will crash GlioVis ...")
                                                    ),
                                                    DT::dataTableOutput(outputId = "knn")
                                             )
                                    ),
                                    
                                    tabPanel(title = "ssGSEA", id = "GSVA", 
                                             p(class = "lead","Single sample Gene Set Enrichment Analysis"),
                                             p(strong("Output format:"), "Subtype calls with the enrichment score."),
                                             br(),
                                             column(width = 9,
                                                    div(class = "busy",  
                                                        p("Calculating, please wait"),
                                                        img(src="Rotating_brain.gif"),
                                                        hr(),
                                                        p("Be patients, switching to another tab will crash GlioVis ...")
                                                    ),
                                                    DT::dataTableOutput(outputId = "gsva")
                                             )
                                    ),
                                    
                                    tabPanel(title = "3-Way", id = "3-Way", 
                                             p(class = "lead","Generate and compare subtype calls by SVM, K-NN and ssGSEA"),
                                             br(),
                                             column(width = 9,
                                                    div(class = "busy",  
                                                        p("Calculating, please wait"),
                                                        img(src="Rotating_brain.gif"),
                                                        hr(),
                                                        p("Be patients, switching to another tab will crash GlioVis ...")
                                                    ),
                                                    DT::dataTableOutput(outputId = "sub3")
                                             )
                                    )
                        )
                      )
                    )
           ),
           
           tabPanel(title ="EstimateME",
                    sidebarLayout(
                      sidebarPanel(width = 3,
                                   helpText("Upload a .csv file with samples in rows and genes expression in columns.",
                                            "The first column should contain the sample ID and should be named 'Sample'"),
                                   radioButtons(inputId = "platformEst", label = strong("Select platform:"), 
                                                choices = c("Affimetrix" = "affymetrix", "Agilent" = "agilent", "Illumina" = "illumina")),
                                   fileInput(inputId = "upEstFile", label = "Choose CSV File",
                                             accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                   hr(),
                                   checkboxInput(inputId = "headerEst", label = "Header", value = TRUE),
                                   radioButtons(inputId = "sepEst", label = "Separator",  
                                                choices = c(Comma = ",", Semicolon = ";", Tab="\t"), selected = ","),
                                   radioButtons(inputId = "quoteEst", label = "Quote", 
                                                choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"'),
                                   hr(),
                                   br(),
                                   conditionalPanel(
                                     condition = "output.finishedEstUploading",
                                     actionButton(inputId = "goEst", label = "Submit", class= "btn-success"),
                                     br()
                                   )
                      ),
                      mainPanel(
                        div(class = "wrap",
                            div(class = "left", p(class = "lead","Estimate of STromal and Immune cells in MAlignant Tumor tissues")),
                            div(class = "right", helpModal(modal_title ="Estimate", link = "helpEst", help_file = includeMarkdown("tools/help_estimate.Rmd")))
                        ),
                        div(class = "busy",  
                            p("Calculating, please wait"),
                            img(src="Rotating_brain.gif") 
                        ),
                        column(width = 7,
                               DT::dataTableOutput(outputId = "estScore")
                        ),
                        column(width = 5,
                               plotOutput(outputId = "purityPlot"),
                               conditionalPanel (condition = "output.purityPlot",
                                                 helpText("In the above plot, the black dot and the gray dash-lines show ESTIMATE tumor purity 
                                                          and 95% prediction interval, respectively. The grey dots represent the background 
                                                          distribution based on 995 samples from TCGA Affymetrix data set.")
                               )
                        )
                      )
                    )
           )
         )
)
