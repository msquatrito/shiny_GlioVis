#' ui.R for Gliovis

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

shinyUI(  
  
  navbarPage(title = strong("GlioVis"), windowTitle = "GlioVis - Visualization Tools for Glioma Datasets", 
             fluid = TRUE, footer = includeHTML("tools/footer.html"), id = "nav",

             source("tools/homeTab.R", local = TRUE)$value,
             source("tools/exploreTab.R", local = TRUE)$value,
             source("tools/toolsTab.R", local = TRUE)$value,
             source("tools/aboutTab.R", local = TRUE)$value,
             source("tools/newsTab.R", local = TRUE)$value
             
  )
)
