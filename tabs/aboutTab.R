# UI-elements for About tab
tabPanel(title = "About", icon = icon("info-circle"),
         div(class="about",
             includeMarkdown("tools/about.Rmd")
         )
)