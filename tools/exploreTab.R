# UI-elements for Explore tab
tabPanel(title = "Explore", icon = icon("picture-o"), id = "explore",
         
         # link to the style.css file. Currently need only for the positioning of the progress bar
         tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
         
         sidebarLayout(
           sidebarPanel(width = 3,
                        img(src = "GlioVis_logo_gray.jpg"),
                        br(),
                        br(),
                        selectInput(inputId = "dataset", label = h4("Dataset"),
                                    choices = datasets,
                                    selected = "TCGA GBM", selectize = TRUE),
                        br(),
                        conditionalPanel(
                          condition = "input.tabCorr == '2genes'",
                          selectizeInput(inputId = "gene", label = h4("Gene"), choices = NULL, selected = NULL, 
                                         options = list(placeholder = "Enter gene, eg: EGFR", plugins = list('restore_on_backspace')))),
                        conditionalPanel(
                          condition = "input.tab1 == 1",
                          br(),
                          selectInput(inputId = "plotTypeSel", label = h4("Plot type"), choices = "", selectize = TRUE),
                          checkboxInput(inputId = "primary", label = "Primary samples only", value = FALSE),
                          br(),
                          helpPopup(title = "Help me pleaseeeeee", content = includeMarkdown("tools/help.Rmd"), 
                                    placement = "bottom", trigger = "click"),
                          hr(),
                          h5("Plot options:"),
                          checkboxInput(inputId = "scale", label = "Scale y axis", value = FALSE),
                          checkboxInput(inputId = "bw", label = "White background", value = FALSE),
                          checkboxInput(inputId = "colBox", label = "Color box", value = FALSE),
                          checkboxInput(inputId = "colStrip", label = "Color stripchart", value = FALSE),
                          conditionalPanel(
                            condition = "input.colStrip == true",
                            uiOutput("colorPoints")
                          ),
                          hr(),
                          h5("Statistic:"),
                          checkboxInput(inputId = "statSummary", label = "Summary statistics", value = FALSE),
                          checkboxInput(inputId = "tukeyHSD", label = "Tukey's HSD", value = FALSE),
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
                          selectInput(inputId = "histologySurv", label = h4("Histology:"), choices = ""),
                          conditionalPanel(
                            condition = "input.histologySurv == 'GBM'",
                            checkboxInput(inputId = "gcimpSurv", label = "Exclude G-CIMP samples", value = FALSE),
                            checkboxInput(inputId = "primarySurv", label = "Exclude Recurrent samples", value = FALSE),
                            br()
                          ),
                          conditionalPanel(
                            condition = "input.histologySurv == 'GBM'",
                            selectInput(inputId = "subtypeSurv", label = h4("Subtype:"), 
                                        choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP"))
                          ),
                          hr(),
                          conditionalPanel(
                            condition = "input.tab1 == 2 & input.tabSurv == 'km'",
                            selectInput(inputId = "cutoff", label = h5("Select cutoff:"), 
                                        choices = c("Use a specific mRNA value", "median", "lower quartile", "upper quartile", "quartiles"),
                                        selected = "median"),
                            conditionalPanel(
                              condition = "input.cutoff == 'Use a specific mRNA value'",
                              br(),
                              uiOutput("numericCutoff"),
                              plotOutput(outputId = "boxmRNA", width = "100%", height = 50),
                              helpText("mRNA expression (log2)")
                            ),
                            hr()
                          ),
                          conditionalPanel(
                            condition = "input.tab1 == 2 & input.tabSurv == 'km'",
                            downloadButton(outputId = "downloadsurvPlot", label = "Download")
                          ),
                          conditionalPanel(
                            condition = "input.tab1 == 2 & input.tabSurv == 'hr'",
                            downloadButton(outputId = "downloadkmPlot", label = "Download")
                          )
                        ),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == '2genes'",
                          br(),
                          selectizeInput(inputId = "gene2", label = h4("Gene 2"), choices ="", 
                                         options = list(placeholder = "Enter gene 2, eg: SOCS2", plugins = list('restore_on_backspace')))),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == 'geneslist'",
                          selectizeInput(inputId = "genelist", label = h4("Genes list"), choices ="", multiple = TRUE,
                                         options = list(placeholder = "Enter genes", plugins = list('remove_button')))),
                        conditionalPanel(
                          condition = "input.tab1 == 3",  
                          br(),
                          selectInput(inputId = "histologyCorr", label = h4("Histology:"), choices = ""),
                          conditionalPanel(
                            condition = "input.histologyCorr == 'GBM'",
                            br(),
                            selectInput(inputId = "subtype", label = h4("Subtype:"), 
                                        choices = c("All", "Classical", "Mesenchymal", "Neural", "Proneural","G-CIMP")))),
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == '2genes'",
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
                        conditionalPanel(
                          condition = "input.tab1 == 3 & input.tabCorr == 'geneslist'",
                        hr(),
                        downloadButton(outputId = "downloadpairsPlot", label = "Download plot")
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
                                    strong("Summary statistics"),
                                    br(),
                                    tableOutput(outputId = "summary"),
                                    hr()
                                  ),
                                  conditionalPanel(
                                    condition = "input.tukeyHSD",
                                    strong("Tukey's Honest Significant Difference (HSD)"),
                                    helpText("The table shows the difference between pairs, the 95% confidence interval and the p-value of the pairwise comparisons:"),
                                    checkboxInput(inputId = "tukeyPlot", label = "Show the results in the plot", value = FALSE),
                                    tableOutput(outputId = "tukeyTest"),
                                    hr()
                                  ),
                                  conditionalPanel(
                                    condition = "input.tTest",
                                    strong("Pairwise t tests"),
                                    helpText("Pairwise comparisons between group levels with corrections for multiple testing (p-values with Bonferroni correction):"),
                                    tableOutput(outputId = "pairwiseTtest")
                                  )
                         ),
                         
                         tabPanel(title = "Survival", icon = icon("user-md"), value = 2,
                                  tabsetPanel(id = "tabSurv",
                                              tabPanel(title = "Km plot",  value = "km",
                                                       plotOutput(outputId = "survPlot", width = 500 , height = 400)
                                              ),                                              
                                              tabPanel(title = "HR plot", value = "hr",
                                                       column(width = 9,
                                                              wellPanel(
                                                                helpText(HTML("<b>IMPORTANT: </b> Currently active only for GBM samples.")),
                                                                helpText(HTML("<b>Note: </b> This is an interactive plot, click on a specific mRNA expression value 
                                                                              to update the survival plot. The blue line represents the current selection.")),
                                                                br(),
                                                                checkboxInput(inputId = "quantile", label = "Show quantiles", value = TRUE),
                                                                # This busy gif does not activate when switching tab
                                                                div(class = "busy",  
                                                                    p("Calculating, please wait"),
                                                                    img(src="ajax-loader.gif") 
                                                                ),
                                                                plotOutput(outputId = "hazardPlot", clickId = "hrClick", width = 500 , height = 400)
                                                              ),
                                                              plotOutput(outputId = "kmPlot", width = 500 , height = 400)
                                                       )
                                                       
                                              )
                                  )
                         ),
                         
                         tabPanel(title = "Correlations", icon = icon("list-alt"),  value = 3,
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
#                                                        ,
#                                                        tableOutput(outputId = "pairsData")
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
                                               uiOutput(outputId = "survPlots", inline = TRUE), # inline = TRUE not working
                                               uiOutput(outputId = "piePlots", inline = TRUE) 
                                             )
                                    )
                                  )
                         )                  
             )
           )
         )
)