# UI-elements for Tools tab
tabPanel(title = "Tools", icon = icon("gear"),
         
         sidebarLayout(fluid = FALSE,
           sidebarPanel(
                        img(src = "GlioVis_tools.jpg", class="responsive-image"),
                        br(),
                        br(),
                        wellPanel(
                          helpText("Upload a .csv file with samples in rows and genes expression in columns.",
                                   "The first column should contain the sample ID and should be named 'Sample'"),
                          conditionalPanel(
                            condition = "input.tabTools == 'SubtypeME'",
                            radioButtons(inputId = "tumorType", label = strong("Training set:"), 
                                         choices = c("GBM" = "gbm", "LGG" = "lgg"))
                          ),
                          conditionalPanel(
                            condition = "input.tabTools == 'EstimateME'",
                            radioButtons(inputId = "platformEst", label = strong("Select platform:"), 
                                         choices = c("Affimetrix" = "affymetrix", "Agilent" = "agilent", "Illumina" = "illumina"))
                          ),
                          conditionalPanel(
                            condition = "input.tabTools == 'DeconvoluteME'",
                            radioButtons(inputId = "platformDec", label = strong("Select platform:"), 
                                         choices = c("Microarray" = "microarray", "RNAseq" = "rnaseq"),inline = T),
                            br(),
                            radioButtons(inputId = "geneListDec", label = strong("Select genelist:"), 
                                         choices = c("Engler et al." = "engler", "Bindea et al." = "galon", 
                                                     "Engler_Bindea" = "galon_engler"),inline = T),
                            br()
                          ),
                          fileInput(inputId = "upFile", label = "Choose Expression data File",
                                    accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                                    ),
                          conditionalPanel(
                            condition = "input.tabTools == 'DeconvoluteME'",
                            checkboxInput(inputId = "deconvPData", label = "Include phenotypic data (eg.: Tumor/Normal)", value = FALSE)
                          ),
                          conditionalPanel(
                            condition = "input.deconvPData",
                            helpText("Upload a .csv file with samples in rows and pData in columns.",
                                     "The first column should contain the sample ID and should be named 'Sample'"),
                            fileInput(inputId = "pDataFile", label = "Choose pData File",
                                      accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                                      ),
                            uiOutput(outputId = "pDataDec")
                          )
                        ),
                        wellPanel(
                          checkboxInput(inputId = "header", label = "Header", value = TRUE)
                        ),
                        wellPanel(
                          radioButtons(inputId = "sep", label = "Separator",  
                                       choices = c(Comma = ",", Semicolon = ";", Tab="\t"), selected = ",", inline = TRUE)
                        ),
                        wellPanel(
                          radioButtons(inputId = "quote", label = "Quote", 
                                       choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"), selected = '"', inline = TRUE)
                        ),
                        conditionalPanel(
                          condition = "input.tabTools == 'SubtypeME'",
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
                        conditionalPanel(
                          condition = "input.tabTools == 'EstimateME'  & output.finishedUploading",
                          actionButton(inputId = "goEst", label = "Submit ESTIMATE", class= "btn-success"),
                          br()
                        ),
                        conditionalPanel(
                          condition = "input.tabTools == 'DeconvoluteME'  & output.finishedUploading",
                          actionButton(inputId = "goDec", label = "Submit Deconvolute", class= "btn-success"),
                          br()
                        )
           ),
           
           mainPanel(
             tabsetPanel(id = "tabTools",
                         tabPanel(title ="SubtypeME",
                                  div(style = "width: 100%; overflow: hidden;",
                                      div(style = "width: 550px; float: left;",p(class = "lead","Classify tumor samples based on mRNA expression profiles")),
                                      div(style = "margin-left: 550px;", helpModal(modal_title ="SubtypeMe", link = "helpClassify", help_file = includeMarkdown("tools/help/help_classify.Rmd")))
                                  ),
                                  
                                  tabsetPanel(id = "tabClassify",   
                                              
                                              tabPanel(title = "SVM", id = "SVM",
                                                       div(style = "width: 100%; overflow: hidden;",
                                                           div(style = "width: 325px; float: left;", p(class = "lead","Supported vector machine learning")),
                                                           div(style = "margin-left: 325px;", helpModal(modal_title ="Supported vector machine", link = "helpSVM", help_file = includeMarkdown("tools/help/help_svm.Rmd")))
                                                       ),
                                                       column(width = 9,
                                                              busy("Be patients, switching to another tab will crash GlioVis ..."),
                                                              dataTableOutput(outputId = "svm")
                                                       )
                                              ),
                                              
                                              tabPanel(title = "K-NN", id = "K-NN",                                            
                                                       div(style = "width: 100%; overflow: hidden;",
                                                           div(style = "width: 290px; float: left;", p(class = "lead","K-nearest neighbors prediction")),
                                                           div(style = "margin-left: 290px;", helpModal(modal_title ="K-nearest neighbors", link = "helpKNN", help_file = includeMarkdown("tools/help/help_knn.Rmd")))
                                                       ),
                                                       column(width = 8,
                                                              busy("Be patients, switching to another tab will crash GlioVis ..."),
                                                              dataTableOutput(outputId = "knn")
                                                       )
                                              ),
                                              
                                              tabPanel(title = "ssGSEA", id = "GSVA",
                                                       div(style = "width: 100%; overflow: hidden;",
                                                           div(style = "width: 415px; float: left;", p(class = "lead","Single sample Gene Set Enrichment Analysis")),
                                                           div(style = "margin-left: 415px;", helpModal(modal_title ="Single sample Gene Set Enrichment Analysis", link = "helpGSEA", help_file = includeMarkdown("tools/help/help_gsea.Rmd")))
                                                       ),
                                                       column(width = 9,
                                                              busy("Be patients, switching to another tab will crash GlioVis ..."),
                                                              dataTableOutput(outputId = "gsva")
                                                       )
                                              ),
                                              
                                              tabPanel(title = "3-Way", id = "3-Way", 
                                                       p(class = "lead","Generate and compare subtype calls by SVM, K-NN and ssGSEA"),
                                                       fluidRow(
                                                         column(width = 10, 
                                                                conditionalPanel(
                                                                  condition = "output.sub3",
                                                                  plotOutput(outputId = "sub3Plot",height = "100%")
                                                                ),
                                                                busy("Be patients, switching to another tab will crash GlioVis ..."),
                                                                dataTableOutput(outputId = "sub3")
                                                         )
                                                       )
                                              )
                                  )
                         ),
                         
                         tabPanel(title ="EstimateME",
                                  div(style = "width: 100%; overflow: hidden;",
                                      div(class = "left", p(class = "lead","Estimate of STromal and Immune cells in MAlignant Tumor tissues")),
                                      div(class = "right", helpModal(modal_title ="Estimate", link = "helpEst", help_file = includeMarkdown("tools/help/help_estimate.Rmd")))
                                  ),
                                  busy(),
                                  splitLayout(cellWidths = c("65%", "35%"),
                                              dataTableOutput(outputId = "estScore"),
                                              plotOutput(outputId = "purityPlot")
                                  )
                         ),
                         
                         tabPanel(title ="DeconvoluteME",
                                  # I've use a slight different HTML style for this Title, I was not able to use the strategy used for the other titles
                                  div(style = "display: inline-block;", p(class = "lead","Deconvolute gene expression profiles into cell-type-specific subprofiles ")),
                                  div(style = "display: inline-block;", helpModal(modal_title ="Deconvolute", link = "helpDeconv", help_file = includeMarkdown("tools/help/help_deconv.Rmd"))),
                                  
                                  tabsetPanel(
                                    tabPanel(title = "Heatmap", id = "deconvHeatmap",
                                             busy(),
                                             plotOutput(outputId = "deconvHeatmap", height = 1000)
                                    ),
                                    tabPanel(title = "Boxplot", id = "deconvBoxPlot",
                                             busy(),
                                             plotOutput(outputId = "deconvBoxPlot", height = "100%")
                                    ),
                                    tabPanel(title = "Scores", id = "deconvScores",
                                             busy(),
                                             dataTableOutput(outputId = "deconvScore")
                                    )
                                    
                                  )
                         )
             )
           )
         )
)

