# UI-elements for About tab
tabPanel(title = "About", icon = icon("info-circle"),
         div(class="about",
             #                           tags$style(type="text/css", '.about {
             #                                         background-image: url("gradient_inv.jpg");
             #                                         background-repeat: no-repeat;
             #                                         background-position: right top;
             #                                         margin-right: 200px;
             #                                       }'),
             includeMarkdown("tools/about.Rmd")
         )
)