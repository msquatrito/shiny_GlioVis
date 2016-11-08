#' install.R for Gliovis

# This file is part of GlioVis
# Copyright (C) Massimo Squatrito
#
# GlioVis is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# 
# GlioVis is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.

#######################################################
# gliovis Installation Instructions: 
# see https://github.com/stan-dev/shinystan/blob/develop/install_shinystan.R

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
  pkg <- c("shiny", "survival", "weights", "googleVis", "tidyverse","shinydashboard", 
           "htmlwidgets","kernlab","devtools","markdown","caret","Cairo", "DT","GGally",
           "Cairo","reshape2", "cgdsr","shinyBS", "gridExtra","survminer")
  
  new.pkg <- pkg[!(pkg %in% installed.packages())]
  
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies=TRUE)
  }
  
  update.packages(pkg[!(pkg %in% new.pkg)])

  # install needed packages from Bioconductor
  bioc <- c("GSVA","limma","DOSE","clusterProfiler")
  if(!(bioc %in% installed.packages())){
    source("http://bioconductor.org/biocLite.R")
    biocLite(bioc) 
  }

  # install needed packages from Github
  devtools::install_github("jokergoo/ComplexHeatmap")
  # devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
  
  message("\n All set. \n You might need to restart R before using GlioVis \n")
  return(invisible(NULL))
}

install_gliovis()
