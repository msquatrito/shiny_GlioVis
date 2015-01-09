library(shiny)
source("global.R")
datasets <- c("TCGA GBM", "TCGA Lgg","Rembrandt", "Gravendeel", "Phillips", "Murat", "Freije", 
              "Reifenberger", "Bao", "Gill", "Gorovets", "Nutt", "Ducray","Grzmil","Donson","Li")

shinyUI(  
  
  navbarPage(title = "GlioVis", windowTitle = "GlioVis - Visualization Tools for Glioma Datasets", fluid = TRUE, responsive = TRUE,
             footer = includeHTML("tools/footer.html"),
             
             source("tools/homeTab.R", local = TRUE)$value,  
             source("tools/exploreTab.R", local = TRUE)$value,
             source("tools/toolsTab.R", local = TRUE)$value,
             source("tools/aboutTab.R", local = TRUE)$value
             
  )
)
