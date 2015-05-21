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
                          selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL,  
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
                          helpModal(modal_title ="Box plots", link = "help", 
                                    help_file = includeMarkdown("tools/help/help.Rmd")),
                          hr(),
                          h4("Plot options:"),
                          div(class="row",
                              div(class="col-xs-6",
                                  numericInput("plot_height", label = "Height (pixels):", min = 100, max = 2000, step = 50, value = 400)
                              ),
                              div(class="col-xs-6",
                                  numericInput("plot_width", label = "Width (pixels):", min = 100, max = 2000, step = 50, value = 500)
                              )
                          ),
                          br(),
                          checkboxInput(inputId = "removeMe", label = "Remove/add a group of samples", value = FALSE),
                          conditionalPanel(
                            condition="input.removeMe",
                            uiOutput(outputId = "removeGp")
                          ),
                          checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                          checkboxInput(inputId = "bw", label = "White background", value = FALSE),
                          checkboxInput(inputId = "colBox", label = "Color box", value = FALSE),
                          checkboxInput(inputId = "colStrip", label = "Color stripchart", value = FALSE),
                          conditionalPanel(
                            condition = "input.colStrip",
                            uiOutput("colorPoints")
                          ),
                          checkboxInput(inputId = "typePoint", label = "Change points appearance", value = FALSE),
                          conditionalPanel(
                            condition="input.typePoint",
                            div(class="row",
                                div(class="col-xs-6",
                                    numericInput(inputId = "point_size",label = "Size", 
                                                 value = 2, min = 0, max = 5, step = 0.5)
                                ),
                                div(class="col-xs-6",
                                    numericInput(inputId = "alpha",label = "Transparency",
                                                 value = 0.5, min = 0, max = 1, step = 0.1)
                                )
                            )
                          ),
                          checkboxInput(inputId = "labelsTitle", label = "Change axis labels", value = FALSE),
                          conditionalPanel(
                            condition="input.labelsTitle",
                            checkboxInput("xaxisLabelAngle", "Rotate sample names", value = FALSE),
                            textInput(inputId = "myXlab", label = "X-axis label:", value = ""),
                            textInput(inputId = "myYlab", label = "Y-axis label:", value = "")
                          ),
                          checkboxInput(inputId = "labelsFont", label = "Change font size", value = FALSE),
                          conditionalPanel(
                            condition="input.labelsFont",
                            div(class="row",
                                div(class="col-xs-6",
                                    numericInput(inputId = "axis_text_size",label = "Axis text",
                                                 value = 12, min = 0, max = 20, step = 1)
                                ),
                                div(class="col-xs-6",
                                    numericInput(inputId = "axis_title_size",label = "Axis title",
                                                 value = 12, min = 0, max = 20, step = 1)
                                )
                            )
                          ),
                          hr(),
                          p(strong("Statistic:")),
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
                            selectInput(inputId = "cutoff", label = "Cutoff:", 
                                        choices = c("Use a specific mRNA value", "median", "lower quartile", "upper quartile", "high vs low", "quartiles"),
                                        selected = "median"),
                            helpModal(modal_title ="Kaplan-Meier", link = "helpKm", help_file = includeMarkdown("tools/help/help_km.Rmd")),
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
                          radioButtons(inputId = "colorBy", label = "Color by:", 
                                       choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                          hr(),
                          radioButtons(inputId = "separateBy", label = "Separate by:", 
                                       choices = c("None" = "none", "Histology" = "Histology", "Subtype" = "Subtype")),
                          hr()
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 5",
                          br(),
                          uiOutput(outputId = "rppaCutoff"),
                          plotOutput(outputId = "boxRppaRNA", width = "100%", height = 50),
                          conditionalPanel(
                            condition = "output.boxRppaRNA",
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
                                              
                                              tabPanel(title = "2-Genes",  value = "2genes",
                                                       tabsetPanel(
                                                         tabPanel(title = "Plot", 
                                                                  plotOutput(outputId = "corrPlot"), 
                                                                  br(),
                                                                  conditionalPanel(
                                                                    condition = "input.statCorr",
                                                                    verbatimTextOutput(outputId = "corrTest")
                                                                  )
                                                         ),
                                                         tabPanel(title = "Data",
                                                                  dataTableOutput(outputId = "corrDataTable", width = 600)
                                                         )
                                                         
                                                       )
                                              ),
                                              
                                              tabPanel(title = "Multiple-Genes", value = "geneslist",
                                                       tabsetPanel(
                                                         tabPanel(title = "Plot",
                                                                  plotOutput(outputId = "pairsPlot")
                                                         ),
                                                         tabPanel(title = "Data",
                                                                  dataTableOutput(outputId = "corrPairsDataTable", width = 600)
                                                         )
                                                       )
                                              )
                                  )
                         ),
                         
                         tabPanel(title = "RPPA", icon = icon("dot-circle-o"), value = 5,
                                  div(style = "width: 100%; overflow: hidden;",
                                      div(style = "width: 550px; float: left;", p(class = "lead","Reverse phase protein array (RPPA) data for TCGA datasets")),
                                      div(style = "margin-left: 550x;", helpModal(modal_title ="RPPA", link = "helpRppa", help_file = includeMarkdown("tools/help/help_rppa.Rmd")))
                                  ),
                                  splitLayout(cellWidths = c("60%", "40%"),
                                              dataTableOutput(outputId = "rppaTable"),
                                              plotOutput(outputId = "rppaPlot", height = 700)
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
                                             dataTableOutput(outputId = "table")
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