# UI-elements for Home tab
tabPanel(title = "Home", icon = icon("home"),
         tagList(
         tags$head(
           includeScript("tools/google-analytics.js"),
           tags$script(type="text/javascript", src = "busy.js"),
           tags$link(rel="shortcut icon", href="./GlioVis_logo_bar.bmp")
           )
         ),
         
         div(id = "home",
                 img(src = "GlioVis_logo_trasnparent.gif", width = 750),
                 hr(),
                 p(class = "lead", "Welcome to", strong("GlioVis"),": a user friendly web based data visualization and analysis application for exploring Glioma datasets."),
                 h4(class = "outer", "How does it work?"),
                 p(class = "outer"," GlioVis is very easy to use:"),
                 tags$ol(
                   tags$li('Select the "Explore" tab'), 
                   tags$li("Choose a dataset"), 
                   tags$li("Enter a Gene Symbol"),
                   tags$li("Select one of the available plots (through the dropdown menu or the specific tab)")
                 ),
                 p(class = "outer", "Available datasets:"),                   
                 source("tools/datsets.table.R",local = TRUE)$value,
                 br(),
                 h4(class = "outer", "Which gene ID can I use?"),
                 p(class = "outer",'Currently only', a("HGNC-approved", href="http://www.genenames.org"), 'protein-coding "Gene Symbols" are supported. Non-coding RNA (miRNA, lncRNA, etc.) are not available.'),
                 h4(class = "outer", "Can I download the plots?"),
                 p(class = "outer",'Yes, all the plots can be downloaded as .pdf files. More file type options are available for the boxplots.'),
                 h4(class = "outer", "Can I download the data?"),
                 p(class = "outer",'Yes, it is', strong("highly recommended"), 'for reproducibility issues. Data can be downloaded at "Explore/Data/Table".'),
                 h4(class = "outer", "What other tools are available?"),
                 p(class = "outer",'SubtypeME: Classify tumor samples based on mRNA expression profiles.'),
                 p(class = "outer",'CorrelateME: Correlate expression of a gene with all the genes in the dataset.'),
                 h4(class = "outer", "Can I use GlioVis results for my publication?"),
                 p(class = "outer", strong("Of course!"), 'If you do so, please include references for the dataset(s) you used and cite:',a("Bowman R. and Squatrito M.", href="#addRef", target="_blank"), " (manuscript in preparation)."),
                 hr(), br(),
                 tags$blockquote(class="pull-right",
                                 tags$p("No great discovery was ever made without a bold guess."), 
                                 tags$small("Isaac Newton")),
                 br(),br(),br(),br(),
                 includeHTML("tools/share.html")
         )
         
)