# UI-elements for Explore tab
tabPanel(title = "Explore", icon = icon("picture-o"), id = "explore",
         
         # link to the style.css file.
         tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
         
         sidebarLayout(
           sidebarPanel(width = 3,
                        img(src = "GlioVis_logo_gray.jpg", class="responsive-image"),
                        br(),
                        br(),
                        selectInput(inputId = "dataset", label = h4("Dataset"),
                                    choices = datasets,
                                    selected = "TCGA GBM", selectize = TRUE),
                        conditionalPanel(
                          condition = "input.tabCorr == '2genes'",
                          selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL, selected = NULL, 
                                         options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace'))),
                          conditionalPanel(
                            condition = "input.tab1 == 1", 
                            radioButtons(inputId ="plotType", label = h4("Plot type"), choices = c("Pre-defined","User-defined"),selected = NULL, inline = T)
                          )
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 1 & input.plotType == 'Pre-defined'",
                          selectInput(inputId = "plotTypeSel", label = NULL, choices = "", selectize = TRUE)
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 1 & input.plotType == 'User-defined'",
                          selectInput(inputId = "plotTypeUserSel", label = NULL, choices = "", selectize = TRUE)
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 1",
                          checkboxInput(inputId = "primary", label = "Primary samples only", value = FALSE),
                          helpModal(modal_title ="Box plots", link = "help", 
                                    help_file = includeMarkdown("tools/help.Rmd")),
#                           helpPopup(title = "Help me pleaseeeeee", content = includeMarkdown("tools/help.Rmd"), 
#                                     placement = "bottom", trigger = "click"),
                          hr(),
                          h5("Plot options:"),
                          checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                          checkboxInput(inputId = "bw", label = "White background", value = FALSE),
                          checkboxInput(inputId = "colBox", label = "Color box", value = FALSE),
                          checkboxInput(inputId = "colStrip", label = "Color stripchart", value = FALSE),
                          conditionalPanel(
                            condition = "input.colStrip",
                            uiOutput("colorPoints")
                          ),
                          hr(),
                          h5("Statistic:"),
                          checkboxInput(inputId = "statSummary", label = "Summary statistics", value = FALSE),
                          checkboxInput(inputId = "tukeyHSD", label = "Tukey's HSD", value = FALSE),
                          checkboxInput(inputId = "tTest", label = "Pairwise t tests", value = FALSE),
                          hr()
                        ),

                        conditionalPanel(
                          condition = "input.tab1 == 2",
                          selectInput(inputId = "histologySurv", label = h4("Histology:"), choices = ""),
                          conditionalPanel(
                            condition = "input.histologySurv == 'GBM'",
                            checkboxInput(inputId = "gcimpSurv", label = "Exclude G-CIMP samples", value = FALSE),
                            checkboxInput(inputId = "primarySurv", label = "Exclude Recurrent samples", value = FALSE)
                          ),
                          conditionalPanel(
                            condition = "input.histologySurv == 'GBM'",
                            selectInput(inputId = "subtypeSurv", label = h4("Subtype:"), 
                                        choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural"))
                          ),
                          conditionalPanel(
                            condition = "input.histologySurv == 'GBM' & input.subtypeSurv == 'All'",
                            checkboxInput(inputId = "allSubSurv", label = "Separate by subtype", value = FALSE)
                          ),
                          hr(),
                          conditionalPanel(
                            condition = "input.tab1 == 2 & input.tabSurv == 'km'",
                            selectInput(inputId = "cutoff", label = h5("Select cutoff:"), 
                                        choices = c("Use a specific mRNA value", "median", "lower quartile", "upper quartile", "high vs low", "quartiles"),
                                        selected = "median"),
                            helpModal(modal_title ="Cutoff selection", link = "helpSurv", 
                                      help_file = includeMarkdown("tools/help_surv.Rmd")),
                            conditionalPanel(
                              condition = "input.cutoff == 'Use a specific mRNA value'",
                              br(),
                              uiOutput("numericCutoff"),
                              plotOutput(outputId = "boxmRNA", width = "100%", height = 50),
                              helpText("mRNA expression (log2). Blue lines represent 25%, 50% and 75% quartiles. Red line represents the current selection.")
                            ),
                            hr()
                          )
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == '2genes'",
                          selectizeInput(inputId = "gene2", label = h4("Gene 2"), choices ="", 
                                         options = list(placeholder = "Enter gene 2, eg: SOCS2", plugins = list('restore_on_backspace')))
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == 'geneslist'",
                          selectizeInput(inputId = "genelist", label = h4("Genes list"), choices ="", multiple = TRUE,
                                         options = list(placeholder = "Enter genes", plugins = list('remove_button')))
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 3",  
                          selectInput(inputId = "histologyCorr", label = h4("Histology:"), choices = ""),
                          conditionalPanel(
                            condition = "input.histologyCorr == 'GBM'",
                            selectInput(inputId = "subtype", label = h4("Subtype:"), 
                                        choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural"))
                          ),
                          hr()
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == '2genes'",
                          h5("Statistic:"),
                          checkboxInput(inputId = "statCorr", label = "Pearson's correlation", value = FALSE),
                          hr(),
                          radioButtons(inputId = "colorBy", label = h5("Color by:"), 
                                       choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                          hr(),
                          radioButtons(inputId = "separateBy", label = h5("Separate by:"), 
                                       choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                          hr()
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 5",
                          br(),
                          uiOutput(outputId = "rppaCutoff"),
                          plotOutput(outputId = "boxRppaRNA", width = "100%", height = 50),
                          conditionalPanel(condition = "output.boxRppaRNA",
                          helpText("mRNA expression (log2). Blue lines represent 25%, 50% and 75% quartiles. Red line represents the current selection."))
                        ),
                        # Allow the user to set the height and width of the plot download.
                        conditionalPanel(
                          condition = "input.tab1 != 4 & input.tab1 != 5",
                          selectInput(inputId = "downloadPlotFileType", label = strong("Select download file type"),
                                      choices = list("PDF"  = "pdf", "BMP"  = "bmp", "JPEG" = "jpeg", "PNG"  = "png")
                          ),
                          strong("Set download image dimensions"),
                          helpText("(units are inches for PDF, pixels for all other formats)"),
                          numericInput(inputId = "downloadPlotHeight", label = "Height (inches)", value = 7, min = 1, max = 100),
                          numericInput(inputId = "downloadPlotWidth", label = "Width (inches)", value = 7, min = 1, max = 100),
                          br(),
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
                            condition = "input.tab1 == 3 & input.tabCorr == '2genes'",
                            downloadButton(outputId = "downloadcorrPlot", label = "Download", class= "btn-primary")
                          ),
                          conditionalPanel(
                            condition = "input.tab1 == 3 & input.tabCorr == 'geneslist'",
                            downloadButton(outputId = "downloadpairsPlot", label = "Download", class= "btn-primary")
                          )
                        )
           ),
           
           
           mainPanel(
             
             h3(textOutput(outputId = "caption")),
             
             tabsetPanel(id = "tab1",
                         
                         tabPanel(title = "Box plots", icon = icon("bar-chart-o"), id = "plots", value = 1,
                                  plotOutput(outputId = "plot"),
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
                         
                         tabPanel(title = "Survival", icon = icon("user-md"), value = 2,
                                  tabsetPanel(id = "tabSurv",
                                              tabPanel(title = "Km plot",  value = "km",
                                                       div(class = "busy",  
                                                           p("Rendering, please wait"),
                                                           img(src="Rotating_brain.gif")),
                                                       plotOutput(outputId = "survPlot", height = "100%")
                                              ),                                              
                                              tabPanel(title = "HR plot", value = "hr",
                                                       column(width = 9,
                                                              wellPanel(
                                                                helpText(HTML('<font color="red"><b>IMPORTANT: </b> Currently active only for GBM samples.</font>')),
                                                                div(class = "busy",  
                                                                    p("Calculating, please wait"),
                                                                    img(src="Rotating_brain.gif") 
                                                                ),
                                                                conditionalPanel(condition = "output.hazardPlot",
                                                                                 helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                                                              to update the survival plot. The blue line represents the current selection.")),
                                                                                 br(),
                                                                                 checkboxInput(inputId = "quantile", label = "Show quantiles", value = TRUE)
                                                                ),
                                                                plotOutput(outputId = "hazardPlot", clickId = "hrClick", width = 500 , height = 400)
                                                              ),
                                                              plotOutput(outputId = "kmPlot", width = 500 , height = 400)
                                                       )
                                                       
                                              )
                                  )
                         ),
                         
                         tabPanel(title = "Correlations", icon = icon("line-chart"),  value = 3,
                                  tabsetPanel(id = "tabCorr",
                                              tabPanel(title = "2-Genes",  value = "2genes",
                                                       plotOutput(outputId = "corrPlot"), 
                                                       br(),
                                                       conditionalPanel(
                                                         condition = "input.statCorr",
                                                         verbatimTextOutput(outputId = "corrTest")
                                                       )
                                              ),
                                              tabPanel(title = "Multiple-Genes", value = "geneslist",
                                                       plotOutput(outputId = "pairsPlot")
                                                       #                ,tableOutput(outputId = "pairsData")
                                              )
                                  )
                         ),
                         
                         tabPanel(title = "RPPA", icon = icon("dot-circle-o"), value = 5,
                                  p(class = "lead","Reverse phase protein array (RPPA) data for TCGA datasets"),
                                  splitLayout(cellWidths = c("60%", "40%"),
                                              DT::dataTableOutput(outputId = "rppaTable"),
                                              plotOutput(outputId = "rppaPlot",height = 700)
                                  )
                         ),
                         
                         tabPanel(title = "Summary", icon =  icon("pie-chart") , value = 4,
                                  tabsetPanel(
                                    
                                    tabPanel(title = "Plots", icon = icon("bar-chart"),
                                             div(class = "busy",  
                                                 p("Rendering report, please wait"),
                                                 img(src="Rotating_brain.gif") 
                                             ),
                                             uiOutput(outputId = "reportPlots")

                                    ),
                                    
                                    tabPanel(title = "Data", icon = icon("table"),
                                             br(),
                                             DT::dataTableOutput(outputId = "table")
                                    ),

                                    tabPanel(title = "Dataset summary plots",
                                             splitLayout(
                                               uiOutput(outputId = "survPlots"),
                                               htmlOutput(outputId = "piePlots") 
                                             )
                                    )
                                  )
                         )                  
             )
           )
         )
)