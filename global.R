#' global.R for Gliovis

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

library(shiny)
library(survival)
# library(weights)
library(gridExtra)
# library(googleVis)
# library(GSVA)
library(GGally)
# library(kernlab)
# library(shinydashboard)
# library(caret)
library(DT)
library(Cairo)
# library(reshape2)
# library(cgdsr)
library(shinyBS)
library(limma)
library(ComplexHeatmap)
library(survminer)
# library(clusterProfiler)
library(tidyverse)
library(plotly)


options(shiny.usecairo=TRUE)

`%then%` <- shiny:::`%OR%`


#######################################
############## Datasets  ##############
#######################################
datasets <- list(
  adult_datasets = c("TCGA_GBM", "TCGA_LGG","TCGA_GBMLGG", "Rembrandt", "Gravendeel","Bao", "Kamoun","Ivy_GAP", "LeeY", "Oh","Phillips",
                     "Gill", "Freije", "Murat","Gorovets", "POLA_Network", "Reifenberger", "Joo","Ducray", "Walsh", "Nutt", "Kwom", "Vital",
                     "Grzmil", "Gleize", "Donson", "Li"),
  pediatric_datasets = c("Northcott_2012","Sturm_2016","Bergthold","Griesinger","Gump","Northcott_2011","Pomeroy",
                         "Johnson","Robinson","Witt", "Kool_2014","Henriquez","Hoffman","Kool_2011","Lambert","Paugh", "deBont","Johann",
                         "Zakrzewski","Sturm_2012","Buczkowicz","Mascelli","Schwartzentruber","Bender"),
  no_surv_dataset = c("Bao","Reifenberger","Gill","Li", "Oh","Ivy_GAP", "Kwom","Walsh","Gleize",
                      "Sturm_2016","Henriquez","Bergthold","Buczkowicz","Mascelli","Lambert","Griesinger",
                      "Zakrzewski","Bender", "deBont","Gump","Johnson","Northcott_2012","Northcott_2011", "Kool_2011",
                      "Johann","Robinson","Kool_2014"),
  rnaseq_datasets = c("TCGA_LGG", "TCGA_GBMLGG", "Bao", "Ivy_GAP","Gill")
)

for (i in unlist(datasets[1:2], use.names = F)) {
  assign(i,readRDS(sprintf("data/datasets/%s.Rds",i))) 
}
rm(i)


#######################################
########## other variables  ###########
#######################################
genes <- readRDS("data/genes.Rds")
lgg.core.samples <- readRDS("data/lgg.core.229samples.Rds")
gene_set_list <- readRDS("data/gene_set_list.Rds")


################################################
##############  Remove NA column  ##############
################################################
.rmNA <- function (df) {
  df <- df[,colSums(is.na(df)) < nrow(df)]
}

##########################################
##############  data table  ##############
##########################################
options(DT.options  = list(lengthMenu = list(c(20, 50, 100, -1), c('20','50','100','All')), 
                           serverSide = FALSE,
                           pagingType = "full",
                           dom = 'lfBrtip',
                           buttons = list('copy', 'print', list(
                             extend = 'collection',
                             buttons = c('csv', 'excel', 'pdf'),
                             text = 'Download'
                           ))
)
)
data_table <- function (df,rownames = FALSE, selection = 'none', filter = "none", options = list()) {
  DT::datatable(df, selection = selection, rownames = rownames, filter = filter, extensions = 'Buttons')
  
}


panel_div <-function(class_type, panel_title, content) {
  HTML(paste0("<div class='panel panel-", class_type,
              "'> <div class='panel-heading'><h3 class='panel-title'>", panel_title,
              "</h3></div><div class='panel-body'>", content,  "</div></div>", sep=""))
}

####################################
##############  busy  ##############
####################################
busy <- function (text = "") {
  div(class = "busy",
      p("Calculating, please wait"),
      img(src="Rotating_brain.gif"),
      hr(),
      p(text)
  )
}

######################################
############## Get HR  ###############
######################################
getHR <- function (df) {
  mRNA <- df[ ,"mRNA"]
  surv.status <- df[ ,"survival_status"]
  surv.time <- df[ ,"survival_month"]
  my.Surv <- Surv(surv.time, surv.status == 1)
  mRNA.values <- mRNA[!is.na(mRNA)]
  # Generate a vector of continuos values, excluding the first an last value
  mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)])
  scan_surv <-function(i, conf.level=95) {
    log.rank <- survdiff(my.Surv ~ mRNA <= i, data = df, rho = 0)
    model <- summary(coxph(my.Surv ~ mRNA <= i))
    HR <- model$conf.int[1]
    HR.lower <- model$conf.int[3]
    HR.upper <- model$conf.int[4]
    data.frame(i, HR, HR.lower, HR.upper, log.rank$obs [1], log.rank$obs [2])
  }
  HRdata <- data.frame (t(sapply(mRNA.values, scan_surv)))
  HRdata <- data.frame (sapply(HRdata,unlist))
  # Exclude groups with less than 10 samples. They don't display properly in the plots(CI too wide)
  HRdata <- subset(HRdata, HRdata[,5] >= 10 & HRdata[,6] >= 10)
}

##########################################
############## hazardPlot  ###############
##########################################
# Function to plot the data from the getHR function
hazardPlot <- function (HRdata, quantile) {
  HRdata[, c("HR", "HR.lower", "HR.upper")] <- log2(HRdata[, c("HR", "HR.lower", "HR.upper")])
  yopt <- median(HRdata[, "HR"])
  ymin <- floor(min(HRdata[, "HR"])) - 1
  if (ymin < yopt-3) ymin <- yopt-3
  if (ymin > -1) ymin <- -2
  ymax <- ceiling(max(HRdata[, "HR"])) + 1
  if (ymax > yopt+3) ymax <- yopt+3
  if (ymax < 1) ymax <- 2
  par(mar = c(5, 4, 4, 2))
  plot(0, 0, type = "n", xlim = range(HRdata[, 1]), ylim = c(ymin - 0.4, ymax + 0.2),
       ylab = "HR with 95% CI", xlab = "mRNA expression (log2)", yaxt="n")
  title(main = "Hazard ratio", line = 3, font.main = 1, cex.main = 1)
  at <- ymin:ymax
  at.pos <- 1:ymax
  at.neg <- abs(ymin):1
  labels <- c(paste(1, 2^at.neg, sep="/"), 1, 2^at.pos)
  axis(2, at = at, labels = labels)
  lines(HRdata[, "i"], HRdata[, "HR"])
  lines(HRdata[, "i"], HRdata[, "HR.lower"], lty=3)
  lines(HRdata[, "i"], HRdata[, "HR.upper"], lty=3)
  points(HRdata[, "i"], rep(ymin - 0.4, length(HRdata[, "i"])), pch="|", col=rgb(0, 0, 0, 0.5))
  abline(h = 0, lty = 2, col = "red")
  if (quantile) {
    q <- quantile(HRdata$i)
    abline(v= q[2:4], col="black")
    axis(3, at = q[2:4], labels = c("25%","50%","75%"),lwd.ticks=0, cex.axis = 0.8, lwd = 0, line = -0.5)
  }
}

###################
## Survival plot ##
###################
get_cutoff <- function(mRNA, cutoff, numeric) {
  mRNA.q <- round(quantile(mRNA, probs=c(0.25, 0.5, 0.75), na.rm = TRUE),2)
  if (cutoff == "quartiles"){
    strat <- cut(mRNA, quantile(mRNA,na.rm = T), include.lowest = TRUE)
  }
  if (cutoff != "quartiles") {
    if (cutoff == "high vs low") {
      strat <- ifelse(mRNA >= mRNA.q[3], "high", ifelse(mRNA <= mRNA.q[1], "low",NA))
    } else {
      cut <- switch(cutoff,
                    "median" = mRNA.q[2],
                    "lower quartile" = mRNA.q[1],
                    "upper quartile" = mRNA.q[3],
                    "Use a specific mRNA value" = numeric)
      f <- function(x) ifelse(x >= cut, c("high"),c("low"))
      strat <- f(mRNA)
    }
  }
  strat
}

survivalPlot <- function (df, gene, group, subtype, cutoff, numeric, censor,risk.table, conf.int, font.legend) {
  # Select a specific Histology
  # if (group != "All") {
  #   df <- filter(df, Histology == group)
  # }
  # Select a specific subtype
  if (subtype != "All") {
    df <- filter(df, Subtype == subtype)
  }
  if(cutoff == "Use a specific mRNA value") {
    main <- paste0("Histology: ", group,
                   "; Subtype: ", subtype,
                   "; Cutoff: ", round(numeric, 2),"\n")
  } else {
    main <- paste0("Histology: ", group,
                   "; Subtype: ", subtype,
                   "; Cutoff: ", cutoff,"\n")
  }
  
  mRNA <- df[ ,"mRNA"]
  surv.status <- df[ ,"survival_status"]
  surv.time <- df[ ,"survival_month"]
  my.Surv <- Surv(time = surv.time, event = surv.status== 1)
  smax <- max(surv.time, na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  df$cutoff_group <- get_cutoff(mRNA, cutoff, numeric)
  expr.surv <- survfit(my.Surv ~ cutoff_group, data=df)
  log.rank <- survdiff(my.Surv ~ cutoff_group, rho = 0, data=df)
  mantle.cox <- survdiff(my.Surv ~ cutoff_group, rho = 1, data=df)
  surv <- data.frame(summary(expr.surv)$table)
  model <- summary(coxph(my.Surv ~ cutoff_group, data=df))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  
  legend.labs = c(sprintf("%s High (n=%s, events=%s, median=%s)", gene, surv$records[1], surv$events[1], surv$median[1]),
                  sprintf("%s Low (n=%s, events=%s, median=%s)", gene, surv$records[2], surv$events[2], surv$median[2]))
  xlegend <- 0.60
  if (cutoff == "quartiles"){
    legend.labs = c(sprintf("1st (n=%s, median=%s)", surv$records[1], surv$median[1]),
                    sprintf("2nd (n=%s, median=%s)", surv$records[2], surv$median[2]),
                    sprintf("3rd (n=%s, median=%s)", surv$records[3], surv$median[3]),
                    sprintf("4th (n=%s, median=%s)", surv$records[4], surv$median[4]))
    xlegend <- 0.70
  }
  p <- survminer::ggsurvplot(data = df, fit = expr.surv, censor = censor, conf.int = conf.int, legend = c(xlegend,0.9), surv.scale = "percent", ylab = "Surviving", xlab = "Survival time (Months)",
                             xlim = c(0,smax),main = main, legend.labs = legend.labs, legend.title = "", font.legend = font.legend, risk.table = risk.table,
                             risk.table.y.text = F, risk.table.y.text.col = T, risk.table.height = 0.4)
  plot <- p$plot
  if (cutoff != "quartiles"){
    plot <- plot + annotate("text", x = xmax, y = c(0.725,0.65,0.575), size = font.legend/3,
                            label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                      sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                      sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p)))
  }
  p$plot <- plot
  print(p)
}

#####################
## Get correlations ##
#####################
# To use to get correlation data (r an p value) on the fly.
# using Hadley suggestion: https://stat.ethz.ch/pipermail/r-help/2008-November/181049.html
getCorr <- function (data, gene, corrMethod) {
  data <- data[ ,9:ncol(data)]
  data <- data[complete.cases(data[,gene]),]
  mRNA <- data[ ,gene, drop = F]
  r <- apply(mRNA, 2, function(x) { apply(data, 2, function(y) { cor(x,y, method = corrMethod) })})
  df <- nrow(mRNA) - 2
  t <- sqrt(df) * r / sqrt(1 - r ^ 2)
  p <- pt(t, df)
  p <- 2 * pmin(p, 1 - p)
  padj <- p.adjust(p, method = "bonferroni")
  corr <- data.frame(row.names(r), round(r,3), signif(p,3),signif(padj,3))
  names(corr) <- c("Gene","r","p.value", "adj.p.value")
  row.names(corr) <- corr$Gene
  corr
}

######################################################
############## 2 genes correlation plot ##############
######################################################
myCorggPlot <- function (df, gene1, gene2, colorBy = "none", separateBy = "none",...) {
  if (separateBy!= "none") {
    df <- df[!is.na(df[,separateBy]), ]
  }
  if (colorBy != "none") {
    df <- df[!is.na(df[,colorBy]), ]
  }
  # empy plot to used in grid.arrange
  empty <- ggplot() + geom_point(aes(1,1), colour="white") + theme_void()
  # scatterplot of x and y variables
  scatter <- ggplot(df,mapping = aes_string(x = gene1, y = gene2)) + theme_bw()
  # marginal density of x - plot on top
  plot_top <- ggplot(df, mapping = aes_string(x = gene1)) + theme_bw() +
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x = element_blank())
  # marginal density of y - plot on the right
  plot_right <- ggplot(df, mapping = aes_string(x = gene2)) + coord_flip() + theme_bw() +
    theme(legend.position = "none",axis.title.y = element_blank(), axis.ticks=element_blank(),axis.text.y=element_blank())
  
  if (colorBy != "none") {
    col <- aes_string(color = colorBy)
    scatter <- scatter + geom_point(col, alpha=0.5) + geom_smooth(col, method = "lm", se = TRUE) + geom_rug(col, alpha = 0.1)
    # grab the legend
    g1 <- ggplotGrob(scatter)
    id.legend <- grep("guide", g1$layout$name)
    legend <- g1[["grobs"]][[id.legend]]
    # scatter without the legend
    scatter <- scatter + theme(legend.position = "none")
    plot_top <- plot_top + geom_density(col, alpha=0.5)
    plot_right <- plot_right + geom_density(col, alpha=0.5)
  }  else {
    scatter <- scatter + geom_point(alpha=0.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1)
    legend <- empty
    plot_top <- plot_top + geom_density(alpha=0.5)
    plot_right <- plot_right + geom_density(alpha=0.5)
  }
  
  if (separateBy != "none") {
    scatter <- scatter + theme(legend.position = "none") + facet_wrap(as.formula(paste("~", separateBy)))
  }  else {
    scatter <- scatter
  }
  
  if (separateBy == "none") {
    # arrange the plots together, with appropriate height and width for each row and column
    grid.draw(grid.arrange(plot_top, legend, scatter, plot_right, ncol=2, nrow=2, widths=c(3, 1), heights=c(1.5, 3)))
  } else {
    print(scatter)
  }
}

######################################
############## kmPlot  ###############
######################################
# Use to plot survival curves, getting the cutoff from the interactive HR plot
kmPlot <- function (cutoff,surv,censor,conf.int,risk.table,font.legend){
  sFit <- survfit(surv)
  sTable <- data.frame(summary(sFit)$table)
  sDiff.log <- survdiff(surv)
  sDiff.mcox <- survdiff(surv,rho = 1)
  log.rank.p <- round(1 - pchisq(sDiff.log$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(sDiff.mcox$chi, df = 1), 4)
  star.log <- weights::starmaker(log.rank.p)
  star.mcox <- weights::starmaker(mantle.cox.p)
  model <- summary(coxph(surv))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  smax <- max(sFit$time,na.rm=TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (90*tmax)/100
  p <- suppressWarnings(
    survminer::ggsurvplot(fit = sFit, censor = censor, conf.int = conf.int, legend = c(0.65,0.9), surv.scale = "percent", ylab = "Surviving", xlab = "Survival time (Months)",
                          xlim = c(0,smax), main = "Kaplan Meier Survival Estimates", legend.title = "", font.legend = font.legend, 
                          legend.labs = c(sprintf("High expr. (n=%s, events=%s, median=%s)", sTable$records[1], sTable$events[1], sTable$median[1]),
                                          sprintf("Low expr. (n=%s, events=%s, median=%s)", sTable$records[2], sTable$events[2], sTable$median[2])),
                          risk.table = risk.table, risk.table.y.text = F, risk.table.y.text.col = T, risk.table.height = 0.4)
  )
  plot <- p$plot
  plot <- plot + annotate("text", x = xmax, y = c(0.725,0.65,0.575), size = font.legend/3,
                          label = c(sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper),
                                    sprintf("%s Log-rank p value= %s", star.log, log.rank.p),
                                    sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p)))
  p$plot <- plot
  print(p)
}

######################################
############## alter_fun_list #######
######################################

alter_fun_list = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue3", col = NA))
  },
  HETLOSS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "cadetblue1", col = NA))
  },
  GAIN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "pink", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

############################################################################
############## Help popup (https://gist.github.com/jcheng5/5913297)  #######
############################################################################
## https://groups.google.com/forum/#!searchin/shiny-discuss/helpPopup/shiny-discuss/ZAkBsL5QwB4/vnmbT47uY7gJ

helpModal <- function(modal_title, link, help_file) {
  sprintf("<div class='modal fade' id='%s' tabindex='-1' role='dialog' aria-labelledby='%s_label' aria-hidden='true'>
          <div class='modal-dialog'>
          <div class='modal-content'>
          <div class='modal-header'>
          <button type='button' class='close' data-dismiss='modal' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
          <h4 class='modal-title' id='%s_label'>%s</h4>
          </div>
          <div class='modal-body'>%s</div>
          </div>
          </div>
          </div>
          <i title='Help' class='fa fa-question-circle' data-toggle='modal' data-target='#%s'></i>",
          link, link, link, modal_title, help_file, link) %>%
    enc2utf8 %>% HTML
}

helpPopup <- function(content, title = NULL) {
  a(href = "#",
    class = "popover-link",
    `data-toggle` = "popover",
    `data-title` = title,
    `data-content` = content,
    `data-html` = "true",
    `data-trigger` = "hover",
    icon("question-circle")
  )
}


#####################################
############## ggally_cor  ##########
#####################################
# To better render ggpairs https://github.com/tonytonov/ggally/blob/master/R/gg-plots.r
ggally_cor <- function(data, mapping, corAlignPercent = 0.6, ...){
  
  xCol <- as.character(mapping$x)
  yCol <- as.character(mapping$y)
  colorCol <- as.character(mapping$colour)
  
  if(length(colorCol) > 0) {
    if(colorCol %in% colnames(data)) {
      rows <- complete.cases(data[,c(xCol,yCol,colorCol)])
    } else {
      rows <- complete.cases(data[,c(xCol,yCol)])
    }
  } else {
    rows <- complete.cases(data[,c(xCol,yCol)])
  }
  
  if(any(!rows)) {
    total <- sum(!rows)
    if (total > 1) {
      warning("Removed ", total, " rows containing missing values")
    } else if (total == 1) {
      warning("Removing 1 row that contained a missing value")
    }
  }
  data <- data[rows, ]
  xVal <- data[,xCol]
  yVal <- data[,yCol]
  
  
  if(length(names(mapping)) > 0){
    for(i in length(names(mapping)):1){
      # find the last value of the aes, such as cyl of as.factor(cyl)
      tmp_map_val <- as.character(mapping[names(mapping)[i]][[1]])
      if(tmp_map_val[length(tmp_map_val)] %in% colnames(data))
        mapping[names(mapping)[i]] <- NULL
      
      if(length(names(mapping)) < 1){
        mapping <- NULL
        break;
      }
    }
  }
  
  final_text <- ""
  if(length(colorCol) < 1)
    colorCol <- "ggally_NO_EXIST"
  
  if(colorCol != "ggally_NO_EXIST" && colorCol %in% colnames(data)) {
    
    txt <- str_c("ddply(data, .(", colorCol, "), summarize, ggally_cor = cor(", xCol,", ", yCol,"))[,c('", colorCol, "', 'ggally_cor')]")
    
    con <- textConnection(txt)
    on.exit(close(con))
    cord <- eval(parse(con))
    
    cord$ggally_cor <- signif(as.numeric(cord$ggally_cor), 3)
    
    # put in correct order
    lev <- levels(data[[colorCol]])
    ord <- rep(-1, nrow(cord))
    for(i in 1:nrow(cord)) {
      for(j in seq_along(lev)){
        if(identical(as.character(cord[i, colorCol]), as.character(lev[j]))) {
          ord[i] <- j
        }
      }
    }
    
    cord <- cord[order(ord[ord >= 0]), ]
    
    cord$label <- str_c(cord[[colorCol]], ": ", cord$ggally_cor)
    
    # calculate variable ranges so the gridlines line up
    xmin <- min(xVal)
    xmax <- max(xVal)
    xrange <- c(xmin-.01*(xmax-xmin),xmax+.01*(xmax-xmin))
    ymin <- min(yVal)
    ymax <- max(yVal)
    yrange <- c(ymin-.01*(ymax-ymin),ymax+.01*(ymax-ymin))
    
    
    p <- ggally_text(
      label   = str_c("Cor : ", signif(cor(xVal,yVal),3)),
      mapping = mapping,
      xP      = 0.5,
      yP      = 0.9,
      xrange  = xrange,
      yrange  = yrange,
      color   = "black",
      ...
    ) +
      #element_bw() +
      theme(legend.position = "none")
    
    xPos <- rep(corAlignPercent, nrow(cord)) * diff(xrange) + min(xrange)
    yPos <- seq(from = 0.9, to = 0.2, length.out = nrow(cord) + 1) * diff(yrange) + min(yrange)
    yPos <- yPos[-1]
    cordf <- data.frame(xPos = xPos, yPos = yPos, labelp = cord$label)
    p <- p + geom_text(
      data=cordf,
      aes(
        x = xPos,
        y = yPos,
        label = labelp,
        color = labelp
      ),
      hjust = 1,
      ...
      
    )
    
    p$type <- "continuous"
    p$subType <- "cor"
    p
  } else {
    # calculate variable ranges so the gridlines line up
    xmin <- min(xVal)
    xmax <- max(xVal)
    xrange <- c(xmin-.01*(xmax-xmin),xmax+.01*(xmax-xmin))
    ymin <- min(yVal)
    ymax <- max(yVal)
    yrange <- c(ymin-.01*(ymax-ymin),ymax+.01*(ymax-ymin))
    
    cor_obj <- cor.test(xVal, yVal)
    cor_signif <- symnum(cor_obj$p.value, corr = FALSE,
                         cutpoints = c(0, .001, .01, .05, .1, 1),
                         symbols = c("***", "**", "*", ".", " "))
    cor_text_size <- 5 + 10*abs(cor_obj$estimate)
    
    p <- ggally_text(
      label = paste(
        signif(cor_obj$estimate, 3),
        cor_signif,
        sep="",collapse=""
      ),
      mapping,
      xP=0.5,
      yP=0.5,
      xrange = xrange,
      yrange = yrange,
      size=cor_text_size,
      ...
    ) +
      #element_bw() +
      theme(legend.position = "none",
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(linetype = "dashed", colour = "black", fill=NA))
    
    p$type <- "continuous"
    p$subType <- "cor"
    p
  }
}

assignInNamespace("ggally_cor", ggally_cor, "GGally")

###############################################
########## majority vote for 3-way sub ########
###############################################
.maj <- function(InVec) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  levels(InVec)[which.max(A)]
}

##################################
########## Estimate score ########
##################################
myEstimateScore <- function (ds, platform = c("affymetrix", "agilent","illumina")) {
  #   stopifnot(is.character(ds) && length(ds) == 1 && nzchar(ds))
  platform <- match.arg(platform)
  data(SI_geneset)
  row.names <- row.names(ds)
  names <- names(ds)
  dataset <- list(ds = ds, row.names = row.names, names = names)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=",
                length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES,
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)
  if (platform != "affymetrix") {
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore", "ImmuneScore",
                              "ESTIMATEScore")
  }
  else {
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884 * x)
    }
    est.new <- NULL
    for (i in 1:length(estimate.score)) {
      est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
      est.new <- rbind(est.new, est_i)
      if (est_i >= 0) {
        next
      }
      else {
        message(paste(sample.names[i], ": out of bounds",
                      sep = ""))
      }
    }
    colnames(est.new) <- c("TumorPurity")
    estimate.t1 <- cbind(estimate.score, est.new)
    x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] <
      0
    estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
    score.data <- rbind(score.data, t(estimate.t1))
    rownames(score.data) <- c("StromalScore", "ImmuneScore",
                              "ESTIMATEScore", "TumorPurity")
  }
  score.data <- data.frame(t(score.data))
}

##################################
########## Plot purity score #####
##################################
.plotPurity <- function (estimate.df, sample, platform = c("affymetrix","agilent", "illumina")) {
  platform <- match.arg(platform)
  convert_row_estimate_score_to_tumor_purity <- function(x) {
    stopifnot(is.numeric(x))
    cos(0.6049872018 + 0.0001467884 * x)
  }
  data(PurityDataAffy)
  Affy.model <- PurityDataAffy
  pred.p <- Affy.model[, 5:7]
  est <- estimate.df[sample, "ESTIMATEScore"]
  purity <- estimate.df[sample, "TumorPurity"]
  max.af <- max(Affy.model$ESTIMATEScore)
  min.af <- min(Affy.model$ESTIMATEScore)
  geMin <- est >= min.af
  leMax <- est <= max.af
  withinMinMax <- geMin && leMax
  xlim <- if (!withinMinMax) {
    adjustment <- 500
    if (geMin) {
      from <- min.af
      to <- est + adjustment
    }
    else {
      from <- est - adjustment
      to <- max.af
    }
    c(from, to)
  }
  else {
    NULL
  }
  plot(Affy.model$tumor.purity ~ Affy.model$ESTIMATEScore,
       Affy.model, main = sample, type = "n", xlab = "ESTIMATE score",
       xlim = xlim, ylab = "Tumor purity", ylim = c(0, 1))
  par(new = TRUE)
  points(Affy.model$ESTIMATEScore, Affy.model$tumor.purity,
         cex = 0.75, col = "lightgrey")
  if (withinMinMax) {
    matlines(Affy.model$ESTIMATEScore, pred.p, lty = c(1,
                                                       2, 2), col = "darkgrey")
  }
  else {
    matlines(Affy.model$ESTIMATEScore, pred.p, lty = c(1,
                                                       2, 2), col = "darkgrey")
    par(new = TRUE)
    curve(convert_row_estimate_score_to_tumor_purity,
          from, to, n = 10000, col = "grey", ylim = c(0,
                                                      1), xlab = "", ylab = "")
  }
  points(est, purity, pch = 19, cex = 1.25)
  abline(h = purity, col = "black", lty = 2)
  abline(v = est, col = "black", lty = 2)
}

############################################
########## ssgsea with permutations ########
############################################
run.ssgsea.GBM <- function (data, number_perms) {
  data <- data.frame(Description = row.names(data),data)
  size <- dim(data[-1])
  num_rows <- size[1]
  num_cols <- size[2]
  num_perm <- number_perms
  random_profile <- data[1]
  for (i in 1:num_perm) {
    a <- data[-1][cbind(seq(1:num_rows), sample(1:num_cols, num_rows, replace = T))]
    random_profile <- cbind(random_profile, a)
    if (i%%100 == 0) {
      print(i)
    }
  }
  
  selected.models <- c("Proneural", "Classical", "Mesenchymal")
  random_result <- .OPAM.apply.model.2(input.ds = random_profile, models.dir = "www/ssgsea_temp/data/", models = selected.models)
  original_result <- .OPAM.apply.model.2(input.ds = data, models.dir = "www/ssgsea_temp/data/", models = selected.models)
  random_result <- random_result[-1]
  random_result <- t(random_result)
  original_result <- original_result[-1]
  original_result <- t(original_result)
  p_result <- original_result
  for (i in 1:dim(original_result)[1]) {
    p_result[i, ] <- colSums(sweep(random_result, 2, original_result[i,]) >= 0)
  }
  
  colnames(p_result) = paste(colnames(p_result), "pval", sep = "_")
  p_result = t(apply(p_result, 1, function(x) {
    (x + 1)/(number_perms + 1)
  }))
  
  results <- data.frame(cbind(round(original_result,2), signif(p_result,3)))
  results$gsea.subtype.call <- as.factor(names(results[4:6]) [apply(results[,4:6],1,which.min)])
  results$gsea.subtype.call <- factor(results$gsea.subtype.call,labels = c("Classical","Mesenchymal", "Proneural"))
  results
}

.OPAM.apply.model.2 <- function (input.ds, input.cls = NA, models.dir, models = "ALL") {
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  m <- data.matrix(input.ds)
  gene.names <- gene.descs <- row.names(input.ds)
  sample.names <- names(input.ds)
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  
  class.labels <- rep(1, Ns)
  class.phen <- "UNDEFINED_PHEN"
  class.list <- rep("U", Ns)
  
  if (models[[1]] == "ALL") {
    file.list <- list.files(models.dir)
    models <- file.list[regexpr(pattern = ".mod", file.list) > 1]
    for (k.model in 1:length(models)) {
      temp <- strsplit(models[k.model], ".mod")
      models[k.model] <- temp[[1]]
    }
    models <- unique(models)
  }
  n.models <- length(models)
  score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  norm.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  model.score.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  probability.matrix <- matrix(0, nrow = n.models, ncol = Ns)
  models.descs <- NULL
  for (model.i in 1:n.models) {
    # print(paste(model.i, "File:", models[model.i]))
    m.file <- paste(models.dir, models[model.i], ".mod", sep = "")
    con <- file(m.file, "r")
    file.content <- readLines(con, n = -1)
    close(con)
    gc()
    len <- length(file.content)
    for (i in 1:len) {
      temp <- unlist(strsplit(file.content[[i]], "\t"))
      len.param <- length(temp)
      if (len.param == 2) {
        param.vals <- temp[2]
      }
      else {
        param.vals <- paste(noquote(temp[2:len.param]), collapse = ",")
        param.vals <- paste("c(", param.vals, ")", sep = "")
      }
      assignment.string <- paste(noquote(temp[1]), " <- ", param.vals, sep = "")
      eval(parse(text = assignment.string))
    }
    
    set.seed(12345)
    if (!exists("sample.norm.type"))
      sample.norm.type <- "rank"
    if (sample.norm.type == "rank") {
      for (j in 1:Ns) {
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- 10000 * m/Ng
    }
    else if (sample.norm.type == "log.rank") {
      for (j in 1:Ns) {
        m[, j] <- rank(m[, j], ties.method = "average")
      }
      m <- log(10000 * m/Ng + exp(1))
    }
    else if (sample.norm.type == "log") {
      m[m < 1] <- 1
      m <- log(m + exp(1))
    }
    
    if (exists("msig.up.genes3"))
      msig.up.genes <- msig.up.genes3
    gene.names.int <- intersect(msig.up.genes, gene.names)
    if (length(gene.names.int) < 2) {
      score.matrix[model.i, ] <- norm.score.matrix[model.i,] <- model.score.matrix[model.i, ] <- probability.matrix[model.i,] <- rep(0, Ns)
      models.descs <- c(models.descs, model.description)
      rm(model.creation.date, input.ds, input.cls, input2.ds,
         input2.cls, target.class, target.class2, model.name,
         model.description, sample.norm.type, marker.disc,
         top.markers.up, top.markers.dn, top.markers.up2,
         top.markers.dn2, statistic, weight, random.seed,
         nperm, link.function, c1, msig.cntrl.genes, msig.up.genes,
         msig.dn.genes, msig.up.genes2, msig.dn.genes2,
         msig.up.genes3, msig.dn.genes3, beta.0, beta.1,
         score, score.up, score.dn)
      next
    }
    
    locs <- match(gene.names.int, gene.names, nomatch = 0)
    msig.up.test <- m[locs, ]
    msig.up.genes.test <- gene.names[locs]
    msig.up.descs.test <- gene.descs[locs]
    msig.up.size.test <- length(locs)
    
    OPAM <- .OPAM.Projection(m, gene.names, Ns, Ng, weight,statistic, msig.up.genes.test, nperm = nperm)
    
    score <- OPAM$ES.vector
    score.matrix[model.i, ] <- score
  }
  
  V.GCT <- data.frame(score.matrix)
  names(V.GCT) <- sample.names
  row.names(V.GCT) <- models
  V.GCT
  
}


.OPAM.Projection <- function (data.array, gene.names, n.cols, n.rows, weight = 0,
                              statistic = "Kolmogorov-Smirnov", gene.set, nperm = 200) {
  ES.vector <- vector(length = n.cols)
  NES.vector <- vector(length = n.cols)
  p.val.vector <- vector(length = n.cols)
  correl.vector <- vector(length = n.rows, mode = "numeric")
  phi <- array(0, c(n.cols, nperm))
  for (sample.index in 1:n.cols) {
    gene.list <- order(data.array[, sample.index], decreasing = T)
    gene.set2 <- match(gene.set, gene.names)
    if (weight == 0) {
      correl.vector <- rep(1, n.rows)
    }
    else if (weight > 0) {
      correl.vector <- data.array[gene.list, sample.index]
    }
    GSEA.results <- .GSEA.EnrichmentScore5(gene.list = gene.list,
                                           gene.set = gene.set2, statistic = statistic, alpha = weight,
                                           correl.vector = correl.vector)
    ES.vector[sample.index] <- GSEA.results$ES
    if (nperm == 0) {
      NES.vector[sample.index] <- ES.vector[sample.index]
      p.val.vector[sample.index] <- 1
    }
    else {
      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:n.rows)
        if (weight == 0) {
          correl.vector <- rep(1, n.rows)
        }
        else if (weight > 0) {
          correl.vector <- data.array[reshuffled.gene.labels, sample.index]
        }
        GSEA.results <- .GSEA.EnrichmentScore5(gene.list = reshuffled.gene.labels,
                                               gene.set = gene.set2, statistic = statistic,
                                               alpha = weight, correl.vector = correl.vector)
        phi[sample.index, r] <- GSEA.results$ES
      }
      if (ES.vector[sample.index] >= 0) {
        pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
        if (length(pos.phi) == 0)
          pos.phi <- 0.5
        pos.m <- mean(pos.phi)
        NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
        s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
        p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
      }
      else {
        neg.phi <- phi[sample.index, phi[sample.index, ] < 0]
        if (length(neg.phi) == 0)
          neg.phi <- 0.5
        neg.m <- mean(neg.phi)
        NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
        s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
        p.val.vector[sample.index] <- ifelse(s == 0,
                                             1/nperm, s)
      }
    }
  }
  return(list(ES.vector = ES.vector, NES.vector = NES.vector,
              p.val.vector = p.val.vector))
}

.GSEA.EnrichmentScore5 <- function (gene.list, gene.set, statistic = "Kolmogorov-Smirnov",
                                    alpha = 1, correl.vector = NULL) {
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  orig.correl.vector <- correl.vector
  if (alpha == 0)
    correl.vector <- rep(1, N)
  correl.vector <- abs(correl.vector)^alpha
  sum.correl <- sum(correl.vector[tag.indicator == 1])
  P0 <- no.tag.indicator/Nm
  F0 <- cumsum(P0)
  Pn <- tag.indicator * correl.vector/sum.correl
  Fn <- cumsum(Pn)
  if (statistic == "Kolmogorov-Smirnov") {
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    }
    else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  }
  else if (statistic == "Cramer-von-Mises") {
    RES <- Fn - F0
    X <- RES^2 * P0
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p)/N)
    ES_n <- sqrt(sum(X_n)/N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    }
    else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
  else if (statistic == "Anderson-Darling") {
    RES <- Fn - F0
    F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1,
                                                          N), F0 * (1 - F0))
    X <- RES^2 * P0/F0_factor
    X_p <- X[RES >= 0]
    X_n <- X[RES < 0]
    ES_p <- sqrt(sum(X_p)/N)
    ES_n <- sqrt(sum(X_n)/N)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(X - max(X_p)))
    }
    else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(X - max(X_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
  else if (statistic == "Zhang_A") {
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn *
                      log(Fn/F0))
    Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl -
                                               1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
    Fn_factor <- ifelse(Fn < 1/sum.correl | Fn > (sum.correl -
                                                    1)/sum.correl, rep(1, N), Fn * (1 - Fn))
    G <- (Fact1 + Fact2) * Pn/Fn_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p)/N
    ES_n <- sum(G_n)/N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    }
    else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
  else if (statistic == "Zhang_C") {
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn *
                      log(Fn/F0))
    Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl -
                                               1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
    F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1,
                                                          N), F0 * (1 - F0))
    G <- (Fact1 + Fact2) * P0/F0_factor
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- sum(G_p)/N
    ES_n <- sum(G_n)/N
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - max(G_p)))
    }
    else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - max(G_n)))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
  else if (statistic == "Zhang_K") {
    RES <- Fn - F0
    Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn *
                      log(Fn/F0))
    Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl -
                                               1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
    G <- Fact1 + Fact2
    G_p <- G[RES >= 0]
    G_n <- G[RES < 0]
    ES_p <- max(G_p)
    ES_n <- max(G_n)
    if (ES_p > ES_n) {
      ES <- signif(ES_p, digits = 5)
      arg.ES <- which.min(abs(G - ES_p))
    }
    else {
      ES <- -signif(ES_n, digits = 5)
      arg.ES <- which.min(abs(G - ES_n))
    }
    return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
  }
  else if (statistic == "area.under.RES") {
    RES <- Fn - F0
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      arg.ES <- which.max(RES)
    }
    else {
      arg.ES <- which.min(RES)
    }
    ES <- sum(RES)
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  }
  else if (statistic == "Wilcoxon") {
    library(exactRankTests)
    seq.index <- seq(1, N)
    gene.set.ranks <- seq.index[tag.indicator == 1]
    gene.set.comp.ranks <- seq.index[tag.indicator == 0]
    W <- wilcox.exact(x = gene.set.ranks, y = gene.set.comp.ranks,
                      alternative = "two.sided", mu = 0, paired = FALSE,
                      exact = F, conf.int = T, conf.level = 0.95)
    ES <- log(1/W$p.value)
    return(list(ES = ES, arg.ES = NULL, RES = NULL, indicator = tag.indicator))
  }
}

# ############################
# ##### Optimal Cutpoint #####
# ############################
# # copy from survminer: https://github.com/kassambara/survminer/blob/master/R/surv_cutpoint.R
plot.surv_cutpoint <- function(x, variables = NULL, ggtheme = theme_classic2(), bins = 30, ...)
{
  
  if(!inherits(x, "surv_cutpoint"))
    stop("x must be an object of class surv_cutpoint.")
  
  data <- x$data
  surv_data <- x$data[, 1:2]
  data <- x$data[, -1*c(1:2), drop = FALSE]
  if(is.null(variables)) variables <- colnames(data)
  data <- data[, variables, drop = FALSE]
  cutpoints <- x$cutpoint[variables,"cutpoint"]
  nvar <- length(variables)
  
  p <- list()
  for(variable in variables){
    max_stat <- x[[variable]]
    
    p_data <- data.frame(
      stats = max_stat$stats,
      cuts = max_stat$cuts,
      grps = .dichotomize(max_stat$cuts, max_stat$estimate)
    )
    
    vline_df <- data.frame(x1 = max_stat$estimate, x2 = max_stat$estimate,
                           y1 = 0, y2 = max(p_data$stats))
    cutpoint_label <- paste0("Cutpoint: ", round(max_stat$estimate,2))
    x1 <- y1 <- x2 <- y2 <- NULL
    max_stat_p <- ggplot(data = p_data, mapping=aes_string("cuts", "stats"))+
      geom_point(aes_string(color = "grps"), shape = 19, size = 0.5)+
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                   data = vline_df, linetype = "dashed", size = 0.5)+
      ggplot2::annotate("text", x = max_stat$estimate, y = 0.5,
                        label = cutpoint_label, size = 4)+
      labs(y = "Standardized Log-Rank Statistic",
           x = variable, title = "Maximally Selected Rank Statistics")
    max_stat_p <- .ggpar(max_stat_p, gggtheme = ggtheme, ...)
    
    distribution <- ggpubr::gghistogram(p_data, x = "cuts", fill = "grps",
                                        main = "Distribution", ylab = "Density",
                                        xlab = "", add_density = TRUE, bins = bins)
    distribution <- .ggpar(distribution, gggtheme = ggtheme, ...)
    res <- list(maxstat = max_stat_p, distribution = distribution)
    
    attr(res, "name") <- variable
    attr(res, "cutpoint") <- max_stat$estimate
    res <- structure(res, class = c("list", "plot_surv_cutpoint"))
    
    p[[variable]] <- res
  }
  p
}
# Helper function
# %%%%%%%%%%%%%%%%%%%%%

.ggpar <- function(p, ggtheme = theme_classic2(),...){
  argmt <- list(...)
  p <- ggpubr::ggpar(p, ggtheme = ggtheme,...)
  if(is.null(argmt$font.x)) p <- p + theme(axis.text.x = element_text(face = "plain"))
  if(is.null(argmt$font.y)) p <- p + theme(axis.text.y = element_text(face = "plain"))
  p
}

# Helper functions
#+++++++++++++++++++++++++++++++++
.dichotomize <- function(x, cutpoint, labels = c("low", "high")){
  grps <- x
  grps[x <= cutpoint] = labels[1]
  grps[x > cutpoint] = labels[2]
  
  grps
}

