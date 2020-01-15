# UI-elements for About tab
tabPanel(title = "About", icon = icon("info-circle"),
         
         panel_div(class_type = "primary",panel_title = "Data Processing",content = includeMarkdown("tools/Processing.Rmd")),
         panel_div(class_type = "primary",panel_title = "Classification",content = includeMarkdown("tools/Classification.Rmd")),
         panel_div(class_type = "primary",panel_title = "Funding",
                  content = HTML ('<a href ="http://www.seveballesteros.com"><img src="seve.jpg" height = "135"  width = "405"></a>')),
         panel_div(class_type = "primary",panel_title = "App info",content = includeMarkdown("tools/App_info.Rmd")),
         panel_div(class_type = "primary",panel_title = "Session info",content = htmlOutput("sessionInfo")),
         panel_div(class_type = "primary",panel_title = "License",content = includeMarkdown("tools/License.Rmd"))

)