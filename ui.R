# ui.R for Gliovis
shinyUI(  
  
  navbarPage(title = "GlioVis", windowTitle = "GlioVis - Visualization Tools for Glioma Datasets", 
             fluid = TRUE, responsive = TRUE, footer = includeHTML("tools/footer.html"),
             
             source("tools/homeTab.R", local = TRUE)$value,  
             source("tools/exploreTab.R", local = TRUE)$value,
             source("tools/toolsTab.R", local = TRUE)$value,
             source("tools/aboutTab.R", local = TRUE)$value,
             source("tools/newsTab.R", local = TRUE)$value
             
  )
)
