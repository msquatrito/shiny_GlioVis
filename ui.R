# ui.R for Gliovis
shinyUI(  
  
  navbarPage(title = strong("GlioVis"), windowTitle = "GlioVis - Visualization Tools for Glioma Datasets", 
             fluid = TRUE, footer = includeHTML("tools/footer.html"), id = "nav",
             
             source("tools/homeTab.R", local = TRUE)$value,  
             source("tools/exploreTab.R", local = TRUE)$value,
             source("tools/toolsTab.R", local = TRUE)$value,
             source("tools/aboutTab.R", local = TRUE)$value,
             source("tools/newsTab.R", local = TRUE)$value, # END navbarMenu
             #### QUIT ####
             tabPanel(title = "Quit", value = "quit", icon = icon("power-off"),
                      h2("Thanks for using GlioVis!"),
                      br(),br(),
                      h5("It's safe to close this browser window.")
             )
             
  )
)
