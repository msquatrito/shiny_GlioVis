# UI-elements for Explore tab
tabPanel(title = "Explore", icon = icon("picture-o"), id = "explore",
         
         # link to the style.css file.
         tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
         
         sidebarLayout(fluid = FALSE,
                       sidebarPanel(
                         img(src = "GlioVis_explore.jpg", class="responsive-image"),
                         br(),
                         br(),
                         wellPanel(                      
                           # dataset
                           selectInput(inputId = "dataset", label = h4("Dataset:"),
                                       choices = datasets,
                                       selected = "TCGA GBM", selectize = TRUE),
                           # genes
                           conditionalPanel(
                             condition = "input.tabCorr != 'corrMany' & input.tab1 != 5 & input.tab1 != 8",
                             selectizeInput(inputId = "gene", label = h4("Gene:"), choices = NULL,  
                                            options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace')))
                           ),
                           conditionalPanel(
                             condition = "input.tab1 == 3 & input.tabCorr == 'corrTwo'",
                             selectizeInput(inputId = "gene2", label = h4("Gene 2:"), choices ="", 
                                            options = list(placeholder = "Enter gene 2, eg: SOCS2", plugins = list('restore_on_backspace')))
                           ),
                           conditionalPanel(
                             condition = "input.tab1 == 3 & input.tabCorr == 'corrMany'",
                             selectizeInput(inputId = "corrGene", label = h4("Genes:"), choices ="", multiple = TRUE,
                                            options = list(placeholder = "Enter two or more genes", plugins = list('remove_button','drag_drop'))),
                             helpText("- or paste -"),
                             textInput(inputId = "genelist_corr", label = h4("Genes list:"),  placeholder = "Paste genes list...")
                           ),
                           conditionalPanel(
                             condition = "input.tab1 == 5",
                             selectizeInput(inputId = "mutGene", label = h4("Gene(s):"), choices ="", multiple = TRUE,
                                            options = list(placeholder = "Enter gene(s), eg: TP53, ATRX", plugins = list('remove_button'))),
                             helpText("- or paste -"),
                             textInput(inputId = "genelist_mut", label = h4("Genes list:"), placeholder = "Paste genes list...")
                           ),
                           
                           # histology and subtype
                           conditionalPanel(
                             condition = "input.tab1 == 2 | input.tab1 == 3 | input.tab1 == 7",
                             selectInput(inputId = "histology", label = h4("Histology:"), choices = ""),
                             selectInput(inputId = "subtype", label = h4("Subtype:"), choices = "")
                           ),
                           
                           # Tab Boxplot
                           conditionalPanel(
                             condition = "input.tab1 == 1", 
                             radioButtons(inputId ="plotType", label = h4("Plot type:"), choices = c("Pre-defined","User-defined"),selected = NULL, inline = T),
                             conditionalPanel(
                               condition = "input.plotType == 'Pre-defined'",
                               selectInput(inputId = "plotTypeSel", label = NULL, choices = "", selectize = TRUE)
                             ),
                             conditionalPanel(
                               condition = "input.plotType == 'User-defined'",
                               selectInput(inputId = "plotTypeUserSel", label = NULL, choices = "", selectize = TRUE)
                             ),
                             checkboxInput(inputId = "removeMe", label = "Remove/add a group of samples", value = FALSE),
                             conditionalPanel(
                               condition="input.removeMe",
                               uiOutput(outputId = "removeGp")
                             ),
                             helpModal(modal_title ="Box plots", link = "help", 
                                       help_file = includeMarkdown("tools/help/help.Rmd"))
                           ),
                           
                           # Tab Survival
                           conditionalPanel(
                             condition = "input.tab1 == 2",
                             conditionalPanel(
                               condition = "input.histology == 'GBM' & input.dataset != 'TCGA GBMLGG'",
                               checkboxInput(inputId = "gcimpSurv", label = "Exclude G-CIMP samples", value = FALSE),
                               checkboxInput(inputId = "primarySurv", label = "Exclude Recurrent samples", value = FALSE)
                             ),
                             conditionalPanel(
                               condition = "input.subtype == 'All'",
                               checkboxInput(inputId = "allSubSurv", label = "Separate by subtype", value = FALSE)
                               
                             ),
                             conditionalPanel(
                               condition = "input.tabSurv == 'km'",
                               selectInput(inputId = "cutoff", label = h4("Cutoff:"), 
                                           choices = c("Use a specific mRNA value", "median", "lower quartile", "upper quartile", "high vs low", "quartiles"),
                                           selected = "median"),
                               helpModal(modal_title ="Kaplan-Meier", link = "helpKm", help_file = includeMarkdown("tools/help/help_km.Rmd")),
                               conditionalPanel(
                                 condition = "input.cutoff == 'Use a specific mRNA value'",
                                 br(),
                                 uiOutput("numericCutoff"),
                                 plotOutput(outputId = "boxmRNA", width = "100%", height = 50),
                                 helpText("mRNA expression (log2). Blue lines represent 25%, 50% and 75% quartiles. Red line represents the current selection.")
                               )
                             )
                           ),
                           
                           # Tab Correlation
                           conditionalPanel(
                             condition = "input.tab1 == 3",
                             helpModal(modal_title ="Correlation", link = "helpCorr", help_file = includeMarkdown("tools/help/help_corr.Rmd"))
                           ),
                           
                           # Tab RPPA
                           conditionalPanel(
                             condition = "input.tab1 == 4",
                             uiOutput(outputId = "rppaCutoff"),
                             plotOutput(outputId = "boxRppaRNA", width = "100%", height = 50),
                             conditionalPanel(
                               condition = "output.boxRppaRNA",
                               helpText("mRNA expression (log2). Blue lines represent 25%, 50% and 75% quartiles. Red line represents the current selection."))
                           ),
                           
  
                           # Tab Differential expression
                           conditionalPanel(
                             condition = "input.tab1 == 7",
                             selectInput(inputId = "DEcutoff", label = h4("Cutoff:"), 
                                         choices = c("high vs low","median", "lower quartile", "upper quartile"),
                                         selected = "high vs low"),
                             helpModal(modal_title ="Differential expression", link = "helpDE", help_file = includeMarkdown("tools/help/help_de.Rmd"))
                           )
                         ),
                         
                         
                         # Tab Boxplot plotting options
                         conditionalPanel(
                           condition = "input.tab1 == 1",
                           wellPanel( 
                             h4("Statistic:"),
                             checkboxInput(inputId = "statSummary", label = "Summary statistics", value = FALSE),
                             checkboxInput(inputId = "tukeyHSD", label = "Tukey's HSD", value = FALSE),
                             checkboxInput(inputId = "tTest", label = "Pairwise t tests", value = FALSE)
                           ),
                           wellPanel( 
                             h4("Plot options:"),
                             p(style = "background-color: #8ABBD7; padding-left:10px; border: 1px solid #E3E3E3;",
                               strong("Plot size")),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput("plot_height", label = "Height (pixels):", min = 100, max = 2000, step = 50, value = 400)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput("plot_width", label = "Width (pixels):", min = 100, max = 2000, step = 50, value = 500)
                                 )
                             ),
                             p(""),      
                             p(style = "background-color: #8ABBD7; padding-left:10px; border: 1px solid #E3E3E3;",
                               strong("Points appearance")),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput(inputId = "point_size",label = "Size", 
                                                  value = 2, min = 0, max = 5, step = 0.5)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput(inputId = "alpha",label = "Transparency",
                                                  value = 0.5, min = 0, max = 1, step = 0.1)
                                 )
                             ),
                             div(class="row",
                                 div(class="col-xs-6",
                                     selectInput("colorP", "Color by:", choices = "None")
                                 ),
                                 div(class="col-xs-6",
                                     selectInput("shapeP", "Shape by:", choices = "None")
                                 )
                             ),
                             p(style = "background-color: #8ABBD7; padding-left:10px; border: 1px solid #E3E3E3;",
                               strong("Axis labels")),
                             div(class="row",
                                 div(class="col-xs-6",
                                     textInput(inputId = "myXlab", label = "X-axis label:", value = "")
                                 ),
                                 div(class="col-xs-6",
                                     numericInput(inputId = "xaxisLabelAngle",label = "X-axis angle",
                                                  value = 0, min = 0, max = 90, step = 15)
                                 )
                             ),
                             textInput(inputId = "myYlab", label = "Y-axis label:", value = "mRNA expression (log2) \n"),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput(inputId = "axis_text_size",label = "Axis text (pt)",
                                                  value = 12, min = 0, max = 20, step = 1)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput(inputId = "axis_title_size",label = "Axis title (pt)",
                                                  value = 12, min = 0, max = 20, step = 1)
                                 )
                             ),
                             p(""),                            
                             p(style = "background-color: #8ABBD7; padding-left:10px; border: 1px solid #E3E3E3;",
                               strong("Miscellaneous")),
                             checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                             checkboxInput(inputId = "bw", label = "White background", value = FALSE),
                             checkboxInput(inputId = "colBox", label = "Color box", value = FALSE)
                           )
                         ),
                         
                         # Tab correlation plotting options
                         conditionalPanel(
                           condition = "input.tab1 == 3 & input.tabCorr == 'corrTwo'",
                           wellPanel(
                             h4("Statistic:"),
                             radioButtons(inputId = "statCorr", label = "Correlation method", choices = c("Pearson", "Spearman", "Kendall"), inline = TRUE)
                           ),
                           wellPanel(
                             h4("Plot options:"),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput("cor_plot_height", label = "Height (pixels):", min = 100, max = 2000, step = 50, value = 400)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput("cor_plot_width", label = "Width (pixels):", min = 100, max = 2000, step = 50, value = 600)
                                 )
                             ),
                             br(),
                             radioButtons(inputId = "colorBy", label = "Color by:", 
                                          choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"),inline =T), 
                             br(),
                             radioButtons(inputId = "separateBy", label = "Separate by:", 
                                          choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype"),inline =T)
                           )
                         ),
                         
                         conditionalPanel(
                           condition = "input.tabCorr == 'corrAll'",
                           wellPanel(
                             selectInput(inputId = "cor", label = h5("Correlation:"), choices = c("All", "Positive", "Negative")),
                             conditionalPanel(
                               condition = "input.cor == 'All'",
                               sliderInput("range", "Range:",min = -1, max = 1, value = c(-0.3,0.3),step = 0.1,round = FALSE)
                             ),
                             radioButtons(inputId = "sign", label = h5("Signficance:"), choices = c(0.05, 0.01),inline = TRUE),
                             radioButtons(inputId = "corrMethod",label = h5("Method:"), choices = c("Pearson", "Spearman"),inline = TRUE)
                           )
                         ),
                         
                         # Tab Mut, Oncoprint plotting options
                         conditionalPanel(
                           condition = "input.tab1 == 5",
                           wellPanel( 
                             h4("Oncoprint options:"),
                             checkboxInput(inputId = "add_cna", label = "Show copy number alterations", value = FALSE),
                             checkboxInput(inputId = "hide_cases", label = "Hide unaltered cases", value = FALSE),
                             checkboxInput(inputId = "column_barplot", label = "Show column bar plot", value = FALSE),
                             checkboxInput(inputId = "row_barplot", label = "Show row bar plot", value = TRUE)
                           )
                         ),
                         
                         # Tab Differential expression options
                         conditionalPanel(
                           condition = "input.tab1 == 7",
                           wellPanel(
                             h4("Filter data:"),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput(inputId = "lfc", label = "Log2 fold change:", min = 0.5, max = 4, step = 0.25, value = 1)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput(inputId = "pvalueDE", label = "p value:", min = 0.001, max = 0.1,step = 0.01, value = 0.05)
                                     # sliderInput(inputId = "pvalueDE", label = "p value:", min = 1e-4, max = 1e-1, step = 1e-2, value = 5e-2)
                                     # radioButtons(inputId = "pvalueDE", label = "p value:",choices = c(0.05, 0.01),inline = TRUE)
                                 )
                             )
                           ),
                           conditionalPanel(
                             condition = "input.tabDE == 'heatmap'",
                           wellPanel(
                             checkboxInput(inputId = "pDataHeatmap", label = "Include annotation to heatmap", value = FALSE),
                             conditionalPanel(
                               condition = "input.pDataHeatmap",
                               selectInput(inputId = "colorSideHeatmap", label = h4("Select group:"), choices = "")
                             )
                           )
                           )
                         ),
                         
                         
                         # Allow the user to set the height and width of the plot download.
                         conditionalPanel(
                           condition = "input.tab1 != 4 & input.tab1 != 5 & input.tab1 != 6 & input.tab1 != 7 & input.tab1 != 8 & input.tabCorr != 'corrAll'",
                           wellPanel( 
                             h4("Download:"),
                             selectInput(inputId = "downloadPlotFileType", label = strong("Select download file type"),
                                         choices = list("PDF"  = "pdf", "BMP"  = "bmp", "JPEG" = "jpeg", "PNG"  = "png")
                             ),
                             strong("Set download image dimensions"),
                             helpText("(units are inches for PDF, pixels for all other formats)"),
                             div(class="row",
                                 div(class="col-xs-6",
                                     numericInput(inputId = "downloadPlotHeight", label = "Height (inches)", 
                                                  value = 7, min = 1, max = 100)
                                 ),
                                 div(class="col-xs-6",
                                     numericInput(inputId = "downloadPlotWidth", label = "Width (inches)", 
                                                  value = 7, min = 1, max = 100)
                                 )
                             ),
                             br(),
                             
                             #Download buttons
                             conditionalPanel(
                               condition = "input.tab1 == 1",
                               downloadButton(outputId = "downloadPlot", label = "Download", class= "btn-primary")
                             ),
                             conditionalPanel(
                               condition = "input.tab1 == 2 & input.tabSurv == 'km'",
                               downloadButton(outputId = "downloadsurvPlot", label = "Download", class= "btn-primary")
                             ),
                             conditionalPanel(
                               condition = "input.tab1 == 2 & input.tabSurv == 'hr'",
                               downloadButton(outputId = "downloadkmPlot", label = "Download", class= "btn-primary")
                             ),
                             conditionalPanel(
                               condition = "input.tab1 == 3 & input.tabCorr == 'corrTwo'",
                               downloadButton(outputId = "downloadcorrPlot", label = "Download", class= "btn-primary")
                             ),
                             conditionalPanel(
                               condition = "input.tab1 == 3 & input.tabCorr == 'corrMany'",
                               downloadButton(outputId = "downloadpairsPlot", label = "Download", class= "btn-primary")
                             )
                           )
                         )
                       ),
                       
                       
                       mainPanel(
                         
                         h3(textOutput(outputId = "caption")),
                         
                         
                         tabsetPanel(id = "tab1",  
                                     
                                     tabPanel(title = "Box plots", icon = icon("bar-chart-o"), value = 1,
                                              
                                              tabsetPanel(id = "tabBoxPlot", 
                                                          
                                                          tabPanel(title = "Plot", 
                                                                   plotOutput(outputId = "plot",height = "100%"),
                                                                   br(),
                                                                   conditionalPanel(
                                                                     condition = "input.statSummary",
                                                                     p(style = "background-color: #F5F5F5; padding-left:10px; border: 1px solid #E3E3E3;", 
                                                                       strong("Summary statistics")),
                                                                     tableOutput(outputId = "summary"),
                                                                     hr()
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.tukeyHSD",
                                                                     p(style = "background-color: #F5F5F5; padding-left:10px; border: 1px solid #E3E3E3;",
                                                                       strong("Tukey's Honest Significant Difference (HSD)")),
                                                                     helpText("The table shows the difference between pairs, the 95% confidence interval and the p-value of the pairwise comparisons:"),
                                                                     checkboxInput(inputId = "tukeyPlot", label = "Show the results in the plot", value = FALSE),
                                                                     tableOutput(outputId = "tukeyTest"),
                                                                     hr()
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.tTest",
                                                                     p(style = "background-color: #F5F5F5; padding-left:10px; border: 1px solid #E3E3E3;",
                                                                       strong("Pairwise t tests")),
                                                                     helpText("Pairwise comparisons between group levels with corrections for multiple testing (p-values with Bonferroni correction):"),
                                                                     tableOutput(outputId = "pairwiseTtest")
                                                                   )
                                                          ),
                                                          
                                                          tabPanel(title = "Data",
                                                                   dataTableOutput(outputId = "filterDataTable", width = 500)
                                                          )
                                              )
                                     ),
                                     
                                     tabPanel(title = "Survival", icon = icon("user-md"), value = 2,
                                              
                                              tabsetPanel(id = "tabSurv",
                                                          
                                                          tabPanel(title = "Kaplan-Meier",  value = "km",
                                                                   tabsetPanel(
                                                                     tabPanel(title = "Plot",
                                                                              p(class = "lead","Kaplan-Meier estimator survival analysis"),
                                                                              plotOutput(outputId = "survPlot", height = "100%")
                                                                     ),
                                                                     tabPanel(title = "Data",
                                                                              dataTableOutput(outputId = "survDataTable", width = 800)
                                                                     )
                                                                   )
                                                          ),  
                                                          
                                                          tabPanel(title = "Hazard ratio", value = "hr",
                                                                   tabsetPanel(
                                                                     tabPanel(title = "Plot",
                                                                              div(style = "width: 100%; overflow: hidden;",
                                                                                  div(style = "width: 120px; float: left;",p(class = "lead","Hazard ratio")),
                                                                                  div(style = "margin-left: 120px;", helpModal(modal_title ="Hazard ratio", link = "helpHR", help_file = includeMarkdown("tools/help/help_hr.Rmd")))
                                                                              ),
                                                                              column(width = 9,
                                                                                     wellPanel(
                                                                                       helpText(HTML('<font color="red"><b>IMPORTANT: </b> Currently active only for GBM samples.</font>')),
                                                                                       busy(),
                                                                                       conditionalPanel(condition = "output.hazardPlot",
                                                                                                        helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                                                              to update the survival plot. The blue line represents the current selection.")),
                                                                                                        br(),
                                                                                                        checkboxInput(inputId = "quantile", label = "Show quantiles", value = TRUE)
                                                                                       ),
                                                                                       plotOutput(outputId = "hazardPlot", click = "hrClick", width = 500 , height = 400)
                                                                                     ),
                                                                                     plotOutput(outputId = "kmPlot", width = 500 , height = 400)
                                                                              )
                                                                              
                                                                     ),
                                                                     
                                                                     tabPanel(title = "Data",
                                                                              dataTableOutput(outputId = "hazardDataTable", width = 600)
                                                                     )
                                                                   )
                                                          )
                                              )
                                     ),
                                     
                                     tabPanel(title = "Correlations", icon = icon("line-chart"),  value = 3,
                                              
                                              tabsetPanel(id = "tabCorr",
                                                          
                                                          tabPanel(title = "Corr-Two",  value = "corrTwo",
                                                                   tabsetPanel(
                                                                     tabPanel(title = "Plot", 
                                                                              plotOutput(outputId = "corrPlot", height = "100%"), 
                                                                              br(),
                                                                              conditionalPanel(
                                                                                condition = "output.corrPlot",
                                                                                p(style = "background-color: #F5F5F5; padding-left:10px; border: 1px solid #E3E3E3;",
                                                                                  strong("Test for association/Correlation between paired samples")),
                                                                                tableOutput(outputId = "corrTest")
                                                                              )
                                                                     ),
                                                                     tabPanel(title = "Data",
                                                                              dataTableOutput(outputId = "corrDataTable", width = 600)
                                                                     )
                                                                     
                                                                   )
                                                          ),
                                                          
                                                          tabPanel(title = "Corr-Many", value = "corrMany",
                                                                   tabsetPanel(
                                                                     tabPanel(title = "Plot",
                                                                              textOutput(outputId = "missingGene"),
                                                                              plotOutput(outputId = "pairsPlot", height = "100%")
                                                                     ),
                                                                     tabPanel(title = "Data",
                                                                              dataTableOutput(outputId = "corrPairsDataTable", width = 600)
                                                                     )
                                                                   )
                                                          ),
                                                          
                                                          tabPanel(title = "Corr-All",  value = "corrAll",
                                                                   p(class = "lead","Correlate expression of a gene with all the genes in the dataset"),
                                                                   busy(),
                                                                   splitLayout(cellWidths = c("60%", "40%"),
                                                                               dataTableOutput(outputId = "corrAllTable"),
                                                                               plotOutput(outputId = "corrAllPlot")
                                                                   )
                                                          )
                                              )
                                     ),
                                     
                                     tabPanel(title = "RPPA", icon = icon("dot-circle-o"), value = 4,
                                              div(style = "width: 100%; overflow: hidden;",
                                                  div(style = "width: 550px; float: left;", p(class = "lead","Reverse phase protein array (RPPA) data for TCGA datasets")),
                                                  div(style = "margin-left: 550x;", helpModal(modal_title ="RPPA", link = "helpRppa", help_file = includeMarkdown("tools/help/help_rppa.Rmd")))
                                              ),
                                              splitLayout(cellWidths = c("60%", "40%"),
                                                          dataTableOutput(outputId = "rppaTable"),
                                                          plotOutput(outputId = "rppaPlot", height = 700)
                                              )
                                     ),
                                     
                                     tabPanel(title = "Mutations",  icon =  icon("adjust"), value = 5,
                                              div(style = "width: 100%; overflow: hidden;",
                                                  div(style = "width: 300px; float: left;", p(class = "lead","Mutation data for TCGA datasets")),
                                                  div(style = "margin-left: 300x;", helpModal(modal_title ="Mutations", link = "helpMut", help_file = includeMarkdown("tools/help/help_mut.Rmd")))
                                              ),
                                              bsAlert("mutAlert"),
                                              plotOutput(outputId = "oncoprint", height = "100%"),
                                              br(),hr(),br(),
                                              busy("Retrieving mutation data"),
                                              dataTableOutput(outputId = "mut")
                                              
                                     ),
                                     
                                     tabPanel(title = "Differential expression", icon = icon("barcode"),  value = 7,
                                              
                                              tabsetPanel(id = "tabDE",
                                                          
                                                          tabPanel(title = "DE Heatmap", value = "heatmap",
                                                                   bsAlert("DEAlert"),
                                                                   busy(),
                                                                   plotOutput(outputId = "DEheatmap", height = 1000)      
                                                          ),
                                                          
                                                          tabPanel(title = "DE Data",
                                                                   busy(),
                                                                   dataTableOutput(outputId = "DETable")        
                                                          ),
                                                          
                                                          tabPanel(title = "Gene Ontology", 
                                                                   p(class = "lead","Gene ontology enrichment analysis"),
                                                                   sidebarLayout(fluid = FALSE,
                                                                                 sidebarPanel(
                                                                                   wellPanel(
                                                                                     radioButtons(inputId = "ont", label = h5("Subontology:"), choices = c( "Biological process" = "BP",
                                                                                                                                                            "Molecular function" = "MF", 
                                                                                                                                                            "Cell compartment" = "CC"))),
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "pvalueCutoff", label = h5("p value cutoff:"), min = 0.0, max = 0.5, step = 0.01, value = 0.05),
                                                                                     sliderInput(inputId = "qvalueCutoff", label = h5("q value cutoff:"), min = 0.0, max = 0.5, step = 0.01, value = 0.05)
                                                                                   ),
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "showCategory", label = h5("# of category to be shown"), min = 1, max = 50, step = 1, value = 10)
                                                                                   ),
                                                                                   actionButton(inputId = "goGO", label = "Submit GO analysis", class= "btn-success"),
                                                                                   helpText("Note: you will have to re-perform the analysis everytime one of the input parameters is changed")
                                                                                 ),
                                                                                 mainPanel(
                                                                                   tabsetPanel(id= "tabGO",
                                                                                               tabPanel(title = "Dotplot", value = "goPlot",
                                                                                                        bsAlert("goAlert"),
                                                                                                        busy(),
                                                                                                        plotOutput(outputId = "enrichGOPlot", width = "800px")
                                                                                               ),
                                                                                               tabPanel(title = "Map",
                                                                                                        busy(),
                                                                                                        plotOutput(outputId = "enrichGOMap", height = "800px")
                                                                                               ),
                                                                                               tabPanel(title = "Data",
                                                                                                        busy(),
                                                                                                        dataTableOutput(outputId = "enrichGOTable")        
                                                                                               )
                                                                                   )
                                                                                 )
                                                                   )
                                                          ),
                                                          
                                                          tabPanel(title = "KEGG", 
                                                                   p(class = "lead","KEGG enrichment analysis"),
                                                                   sidebarLayout(fluid = FALSE,
                                                                                 sidebarPanel(
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "pvalueCutoffKegg", label = h5("p value cutoff:"), min = 0.0, max = 0.5, step = 0.01, value = 0.05),
                                                                                     sliderInput(inputId = "qvalueCutoffKegg", label = h5("q value cutoff:"), min = 0.0, max = 0.5, step = 0.01, value = 0.1)
                                                                                   ),
                                                                                   wellPanel(
                                                                                     sliderInput(inputId = "showCategoryKegg", label = h5("# of category to be shown"), min = 1, max = 50, step = 1, value = 10)
                                                                                   ),
                                                                                   actionButton(inputId = "goKegg", label = "Submit KEGG analysis", class= "btn-success"),
                                                                                   helpText("Note: you will have to re-perform the analysis everytime one of the input parameters is changed")
                                                                                 ),
                                                                                 mainPanel(
                                                                                   tabsetPanel(id= "tabKegg",
                                                                                               tabPanel(title = "Dotplot", value = "kggPlot",
                                                                                                        bsAlert("keggAlert"),
                                                                                                        busy(),
                                                                                                        plotOutput(outputId = "enrichKeggPlot", width = "800px")
                                                                                               ),
                                                                                               tabPanel(title = "Map",
                                                                                                        busy(),
                                                                                                        plotOutput(outputId = "enrichKeggMap", height = "800px")
                                                                                               ),
                                                                                               tabPanel(title = "Data",
                                                                                                        busy(),
                                                                                                        dataTableOutput(outputId = "enrichKeggTable")        
                                                                                               )
                                                                                   )
                                                                                 )
                                                                   )
                                                          )
                                              )
                                     ),
                                     
                                     tabPanel(title = "Summary", icon =  icon("pie-chart"), value = 6,
                                              
                                              tabsetPanel(
                                                
                                                tabPanel(title = "Plots", icon = icon("bar-chart"),
                                                         uiOutput(outputId = "reportPlots")
                                                ),
                                                
                                                tabPanel(title = "Data", icon = icon("table"),
                                                         br(),
                                                         dataTableOutput(outputId = "table")
                                                ),
                                                
                                                tabPanel(title = "Dataset summary plots",
                                                         splitLayout(
                                                           uiOutput(outputId = "survPlots"),
                                                           htmlOutput(outputId = "piePlots") 
                                                         )
                                                )
                                              )
                                     ),
                                     
                                     tabPanel(title = "Dataset download", icon =  icon("download"), value = 8,
                                              p(class = "lead","Download complete dataset"),
                                              br(),
                                              downloadButton(outputId = "downloadExpressionData", label = "Download expression data", class= "btn-primary"),
                                              p(""),
                                              downloadButton(outputId = "downloadpData", label = "Download pheno data", class= "btn-primary")
                                     )
                                     
                         )
                       )
         )
)