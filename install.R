#######################################################
# gliovis Installation Instructions: 
# https://github.com/stan-dev/shinystan/blob/develop/install_shinystan.R

# 1) Restart R
# 2) Run install_gliovis(), included below, to install

#######################################################

install_gliovis <- function() {
  
  if (getRversion() < '3.1.1') stop("GlioVis requires R version 3.1.1 or greater.")
  
  msg <- "Note: this will install or update packages needed to run GlioVis Do you want to continue?"
  continue <- select.list(choices = c("Yes", "No"), title = msg, graphics = FALSE)
  if (continue == "No") {
    message("Installation canceled by user.")
    return(invisible(NULL))
  }
  
  # install needed packages from CRAN
  pkg <- c("shiny", "survival", "weights", "googleVis", "dplyr", "ggplot2","shinydashboard", 
           "htmlwidgets","kernlab","devtools","GGally","markdown","caret","Cairo","broom",
           "Cairo","reshape2","heatmap3","RColorBrewer","scales","htmlwidgets","DT")
  
  new.pkg <- pkg[!(pkg %in% installed.packages())]
  
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies=TRUE)
  }
  
  update.packages(pkg[!(pkg %in% new.pkg)])
  
  if (!require("GSVA")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("GSVA")
  }
  
  if (packageVersion('gridExtra') > '0.9.1') {
    packageurl <- "http://cran.r-project.org/src/contrib/Archive/gridExtra/gridExtra_0.9.1.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
  }
  
  # if (!require("estimate")){
  #   library(utils)
  #   mdacc <- local({
  #     rvers <- getRversion()
  #     repos.hostname <- "bioinformatics.mdanderson.org"
  #     sprintf("http://%s/OOMPA/%s",
  #             repos.hostname,
  #             paste(rvers$maj, rvers$min, sep="."))
  #   })
  #   install.packages("estimate", repos=mdacc, dependencies=TRUE)
  # }
  
  # install needed packages from GitHub
#   if (!require("shinydashboard")) devtools::install_github("rstudio/shinydashboard")
#   if (!requireNamespace('htmlwidgets') || packageVersion('htmlwidgets') <= '0.3.2')
#     devtools::install_github('ramnathv/htmlwidgets')
#   if (!require("DT")) devtools::install_github("rstudio/DT")
  
  message("\n All set. \n You might need to restart R before using GlioVis \n")
  return(invisible(NULL))
}
