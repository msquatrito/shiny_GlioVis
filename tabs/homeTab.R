# UI-elements for Home tab
tabPanel(title = "Home", icon = icon("home"),
         # absolutePanel(top = 100, right = 20, width = 400, draggable =TRUE, cursor = "move",
         #               div(style = "border: 1px solid; padding:10px; border-radius: 15px;",
         #                   # span(class="label label-danger","NEWS"),
         #                 p(style = "background-color: #68cbd0; border-radius: 15px; padding-left:10px; border: 1px solid #E3E3E3;",
         #                   strong("NEWS")),
         #                 includeMarkdown("tools/news_home.Rmd")
         #               ), style = "opacity: 0.92"),
         tagList(
           tags$head(
             includeScript("tools/google-analytics.js"),
             tags$script(type="text/javascript", src = "busy.js"),
             tags$script(type="text/javascript", src = "slider.js"),
             tags$link(rel="shortcut icon", href="./GlioVis_logo_bar.bmp"),
             tags$script(type="text/javascript", "var switchTo5x=true"),
             tags$script(type="text/javascript", src="http://w.sharethis.com/button/buttons.js"),
             tags$script(type="text/javascript",'stLight.options({publisher: "675b3562-a081-470a-9fc4-3dd6a712209d", doNotHash: true, doNotCopy: true, hashAddressBar: false})')
           )
         ),

         div(id = "home",
             div(class="pull-right",
                 span(class='st_twitter', displayText='Tweet'),
                 span(class='st_linkedin', displayText='LinkedIn'),
                 span(class='st_facebook', displayText='Facebook'),
                 span(class='st_email', displayText='Email')
             ),
             br(),
             a(href ="https://www.cnio.es/eventos/index.asp?ev=1&cev=136",div(class="simple-ss",id="simple-ss")),
             # img(src = "new_logo_transparent.png", width = 900),
             div(class="intro-divider"),
             p(class = "lead", "Welcome to", strong("GlioVis"),": a user friendly web application for data visualization and analysis to explore brain tumors expression datasets."),
             div(style = "border: 1px solid; padding:10px; border-radius: 15px;",
                 p(style = "background-color: #68cbd0; border-radius: 15px; padding-left:10px; border: 1px solid #E3E3E3;",
                   strong("NEWS")),
                 includeMarkdown("tools/news_home.Rmd")
             ),
             h4(class = "outer", "How does it work?"),
             p(class = "outer"," GlioVis is very easy to use:"),
             tags$ol(
               tags$li('Select the "Explore" tab'),
               tags$li("Choose a dataset"),
               tags$li("Enter a Gene Symbol"),
               tags$li("Select one of the available plots (through the dropdown menu or a specific tab)")
             ),
             h4(class = "outer", "What type of plots are available?"),
             p(class = "outer",'For each dataset it is available a list of common "pre-defined" plots (see the table below) and
                                a list of dataset-specific "user-defined" box plots. An overview of all the box plots for a
                                given dataset can be found at "Explore/Summary/Plots".'),
             p(class = "outer", "Available datasets:"),
             tabsetPanel(id = "datasets", #type = "pills",
                         tabPanel(title=p("Adult"),
                                  source("tools/adult_datasets_table.R",local = TRUE)$value),
                         tabPanel(title= p("Pediatric"),
                                  #title= p("Pediatric",span(class="label label-danger","New")),
                                  #HTML('<p>Pediatric <span class="label label-danger">New</span></p>'),
                                  source("tools/pediatric_datasets_table.R",local = TRUE)$value)
             ),
             br(),
             h4(class = "outer", "Which gene ID can I use?"),
             p(class = "outer",'Currently only', a("HGNC-approved", href="http://www.genenames.org"), 'protein-coding "Gene Symbols" are supported.
               Non-coding RNA (miRNA, lncRNA, etc.) are not available.'),
             h4(class = "outer", "Can I download the plots?"),
             p(class = "outer",'Yes, all the plots can be downloaded in various formats: .pdf, .bmp, .jpeg or .jpg.'),
             h4(class = "outer", "Can I download the data?"),
             p(class = "outer",'Yes, it is', strong("highly recommended"), 'for reproducibility issues. For each plot there is a "Data" tab containing
               the actual data used to generate the plot. The data can also be downloaded at "Explore/Summary/Data".'),
             h4(class = "outer", "What other tools are available?"),
             p(class = "outer",'SubtypeME: Classify tumor samples based on mRNA expression profiles.'),
             p(class = "outer",'EstimateME: Estimate of STromal and Immune cells in MAlignant Tumor samples',
               a("(Yoshihara K. et al., 2013).", href="http://www.ncbi.nlm.nih.gov/pubmed/24113773", target="_blank")),
             p(class = "outer",'DeconvoluteME: Deconvolute gene expression profiles from heterogeneous tissue samples into cell-type-specific subprofiles '),
             h4(class = "outer", "Can I use GlioVis results for my publication?"),
             p(class = "outer", strong("Of course!"), 'If you do so, please include references for the dataset(s) you used and cite:', a("Bowman R. and Squatrito M.", href="#addRef", target="_blank"), " (manuscript in preparation)."),
             p(class = "outer","Please adhere to the",  a("TCGA publication guidelines", href="http://cancergenome.nih.gov/publications/publicationguidelines"),
               "when using TCGA data in your publications."),
             hr(), br(),
             tags$blockquote(class="pull-right",
                             tags$p("No great discovery was ever made without a bold guess."),
                             tags$small("Isaac Newton")),
             br(),br(),br(),br(),
             # includeHTML("tools/share.html"),
             includeHTML("tools/disqus.html")
         )

)
