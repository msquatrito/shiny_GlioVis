pkg <- c("shiny", "survival", "weights", "googleVis", "dplyr", "ggplot2","gridExtra", 
        "htmlwidgets","kernlab","devtools","GGally","markdown","caret","Cairo")
new.pkg <- pkg[!(pkg %in% installed.packages())]
if (length(new.pkg)) {
  install.packages(new.pkg, dependencies=TRUE)
}
if (!require("GSVA")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GSVA")
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

if (!require("shinydashboard")) devtools::install_github("rstudio/shinydashboard")
if (!requireNamespace('htmlwidgets') || packageVersion('htmlwidgets') <= '0.3.2')
  devtools::install_github('ramnathv/htmlwidgets')
if (!require("DT")) devtools::install_github("rstudio/DT")
