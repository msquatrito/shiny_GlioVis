library(shiny)
library(survival)
library(weights)
library(ggplot2)
library(gridExtra)
library(googleVis)
library(dplyr)
library(parallel)
library(GSVA)
library(GGally)
library(kernlab)
library(shinydashboard)
library(caret)
library(DT)
library(Cairo)
options(shiny.usecairo=TRUE)


#######################################
############## Datasets  ##############
#######################################
datasets <- c("TCGA GBM", "TCGA Lgg","Rembrandt", "Gravendeel", "Phillips", "Murat", "Freije", 
              "Reifenberger", "Bao", "Gill", "Gorovets", "Nutt", "Ducray","Grzmil","Donson","Li","Vital",
              "Joo","Oh")

noSurvDataset <- c("Bao","Reifenberger","Gill","Li", "Oh")

gbm.tcga <- readRDS("data/TCGA.GBM.Rds")
lgg.tcga <- readRDS("data/TCGA.LGG.Rds")
rembrandt <- readRDS("data/Rembrandt.Rds")
freije <- readRDS("data/Freije.Rds")
gravendeel <- readRDS("data/Gravendeel.Rds")
murat <- readRDS("data/Murat.Rds")
phillips <- readRDS("data/Phillips.Rds")
reifenberger <- readRDS("data/Reifenberger.Rds")
bao <- readRDS("data/Bao.Rds")
gill <- readRDS("data/Gill.Rds")
gorovets <- readRDS("data/Gorovets.Rds")
nutt <- readRDS("data/Nutt.Rds")
ducray <- readRDS("data/Ducray.Rds")
grzmil <- readRDS("data/Grzmil.Rds")
donson <- readRDS("data/Donson.Rds")
li <- readRDS("data/Li.Rds")
vital <- readRDS("data/Vital.Rds")
joo <- readRDS("data/Joo.Rds")
oh <- readRDS("data/Oh.Rds")

#######################################
########## other variables  ###########
#######################################
genes <- readRDS("data/genes.Rds")
gene_names <- as.character(genes[,"Gene"])
gbm.subtype.list <- readRDS("data/subtype_list.Rds")
gbm.core.samples <- readRDS("data/TCGA.core.345samples.Rds")
lgg.core.samples <- readRDS("data/lgg.core.460samples.Rds")

#######################################
############## plotList  ##############
#######################################
plotList <- list("TCGA GBM" = c("Histology", "Copy_number", "Subtype", "CIMP_status", "Recurrence"),
                 "TCGA Lgg" = c("Histology", "Grade", "Copy_number", "Subtype"),
                 "Rembrandt" = c("Histology", "Grade", "Subtype", "CIMP_status"),
                 "Gravendeel" = c("Histology", "Grade", "Subtype", "CIMP_status"),
                 "Phillips" = c("Histology", "Grade", "Subtype", "Recurrence", "CIMP_status"),
                 "Murat" = c("Histology", "Subtype", "Recurrence", "CIMP_status"),
                 "Freije" = c("Histology", "Grade", "Subtype", "CIMP_status"),
                 "Reifenberger" = c("Subtype", "CIMP_status"),
                 "Bao" = c("Histology", "Subtype", "Recurrence", "CIMP_status"),
                 "Gill" = c("Histology", "Subtype", "CIMP_status"),
                 "Gorovets" = c("Histology", "Grade", "Subtype"),
                 "Nutt" = c("Histology", "Subtype", "CIMP_status"),
                 "Ducray" = c("Subtype", "CIMP_status"),
                 "Grzmil"= c("Histology", "Subtype", "CIMP_status"),
                 "Donson"= c("Histology", "Subtype"),
                 "Li" = c("Subtype", "CIMP_status"),
                 "Vital" = c("Histology", "Grade", "Subtype"),
                 "Joo" = c("Histology", "Subtype", "Recurrence", "CIMP_status"),
                 "Oh" = c("Recurrence", "Subtype", "CIMP_status"))


################################################
##############  Remove NA column  ##############
################################################
rmNA <- function (df) {
  df <- df[,colSums(is.na(df)) < nrow(df)]
}

######################################
############## Get HR  ###############
######################################
getHR <- function (df) {
  mRNA <- df[ ,"mRNA"]
  surv.status <- df[ ,"status"]
  surv.time <- df[ ,"survival"]
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
survivalPlot <- function (df, gene, group, subtype, cutoff, numeric) {
  # Select a specific subtype
  if (group == "GBM" & subtype != "All") {
    df <- filter(df, Subtype == subtype)
  }
  mRNA <- df[ ,"mRNA"]
  surv.status <- df[ ,"status"]
  surv.time <- df[ ,"survival"]
  my.Surv <- Surv(time = surv.time, event = surv.status== 1)
  smax <- max(surv.time, na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  xmax <- (95*tmax)/100
  mRNA.q <- round(quantile(mRNA, probs=c(0.25, 0.5, 0.75), na.rm = TRUE),2)
  
  if(cutoff == "Use a specific mRNA value") {
    main <- paste("Histology: ", group, 
                  "; Subtype: ", subtype,
                  "; Cutoff: ", round(numeric, 2), sep = "") 
  } else {
    main <- paste("Histology: ", group, 
                  "; Subtype: ", subtype,
                  "; Cutoff: ", cutoff, sep = "")
  }
  
  if (cutoff != "quartiles") {
    if (cutoff == "high vs low") {
      strat <- ifelse(mRNA >= mRNA.q [3], "high", ifelse(mRNA <= mRNA.q [1], "low",NA))
    } else {
      cut <- switch(cutoff, 
                    "median" = mRNA.q[2],
                    "lower quartile" = mRNA.q [1],
                    "upper quartile" = mRNA.q [3],
                    "Use a specific mRNA value" = numeric)
      f <- function(x) ifelse(x >= cut, c("high"),c("low"))
      strat <- f(mRNA)
    }
    expr.surv <- survfit(my.Surv ~ strat, conf.type = "none")
    log.rank <- survdiff(my.Surv ~ strat, rho = 0)
    mantle.cox <- survdiff(my.Surv~ strat, rho = 1)
    surv <- data.frame(summary(expr.surv)$table)
    model <- summary(coxph(my.Surv ~ strat))
    HR <- round(model$conf.int[1],2)
    HR.lower <- round(model$conf.int[3],2)
    HR.upper <- round(model$conf.int[4],2)
    log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
    mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
    star.log <- starmaker(log.rank.p)
    star.mcox <- starmaker(mantle.cox.p)
    plot(expr.surv, xlab = "Survival time (Months)", ylab = "% Surviving", yscale = 100, xlim = c(0,smax),
         main = main, col = c("red", "blue"), mark.time = FALSE)
    legend("topright", legend = c(sprintf("%s High, (n=%s, events=%s, median=%s)", gene, surv$records[1], surv$events[1], surv$median[1]), 
                                  sprintf("%s Low, (n=%s, events=%s, median=%s)", gene, surv$records[2], surv$events[2], surv$median[2])),
           col= c("red", "blue"), lty = 1, cex = 1)
    text (xmax, 0.725, sprintf("HR = %s, (%s - %s)",HR, HR.lower, HR.upper), cex = 1)
    text (xmax, 0.65, sprintf("%s Log-rank p value= %s", star.log, log.rank.p), cex = 1)
    text (xmax, 0.575, sprintf("%s Wilcoxon p value= %s",star.mcox, mantle.cox.p), cex = 1)
  }
  
  if (cutoff == "quartiles"){
    exprcat <- cut(mRNA, quantile(mRNA,na.rm = T), data=df)
    expr.surv <- survfit(my.Surv ~ strata(exprcat), data=df, conf.type="none")
    z <- data.frame(summary(expr.surv)$table) 
    plot(expr.surv, xlab="Months", ylab="% Surviving", yscale = 100, xlim = c(0,smax), 
         main = main, col= c(1:4), mark.time=FALSE)
    legend("topright", title = "Quantile", 
           legend = c(sprintf("1st (n=%s, median=%s)", z$records[1], z$median[1]),
                      sprintf("2nd (n=%s, median=%s)", z$records[2], z$median[2]),
                      sprintf("3rd (n=%s, median=%s)", z$records[3], z$median[3]),
                      sprintf("4th (n=%s, median=%s)", z$records[4], z$median[4])),
           col= c(1:4), lty=1, cex=1)
  }
  
}

#####################
## Get correlations ##
#####################
# To use to get correlation data (r an p value) on the fly. 
# using Hadley suggestion: https://stat.ethz.ch/pipermail/r-help/2008-November/181049.html
getCorr <- function (data, gene, histology, corrMethod) {
  if (histology != "All") {
    data <- filter (data, Histology == histology)
  }
  data <- data[ ,9:ncol(data)]
  mRNA <- data[ ,gene, drop = F]
  r <- apply(mRNA, 2, function(x) { apply(data, 2, function(y) { cor(x,y, method = corrMethod) })})
  df <- nrow(mRNA) - 2
  t <- sqrt(df) * r / sqrt(1 - r ^ 2)
  p <- pt(t, df)
  p <- 2 * pmin(p, 1 - p)
  padj <- p.adjust(p, method = "bonferroni")
  corr <- data.frame(row.names(r), round(r,3), round(p,10),round(padj,10))
  names(corr) <- c("Gene","r","p.value", "adj.p.value")
  corr
}

######################################################
############## 2 genes correlation plot ##############
######################################################
myCorggPlot <- function (df, gene1, gene2, histo = "All", subtype = "All", colorBy = "none", separateBy = "none",...) {
  if (histo != "All") {
    df <- filter (df, Histology == histo)
  } 
  if (histo == "GBM" & subtype != "All") {
    df <- filter (df, Subtype == subtype)
  }
  #  empy plot to used in grid.arrange 
  empty <- ggplot() + geom_point(aes(1,1), colour="white") + 
    theme(plot.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.border = element_blank(), 
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank()
    )
  # scatterplot of x and y variables
  scatter <- ggplot(df,mapping = aes_string(x = gene1, y = gene2)) + theme(legend.position=c(1,1),legend.justification=c(1,1)) 
  # marginal density of x - plot on top
  plot_top <- ggplot(df, mapping = aes_string(x = gene1)) + 
    theme(legend.position = "none", axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x = element_blank())
  # marginal density of y - plot on the right
  plot_right <- ggplot(df, mapping = aes_string(x = gene2)) + coord_flip() + 
    theme(legend.position = "none",axis.title.y = element_blank(), axis.ticks=element_blank(),axis.text.y=element_blank())
  
  if (colorBy != "none") {
    col <- aes_string(color = colorBy)
    scatter <- scatter + geom_point(col, alpha=.5) + geom_smooth(col, method = "lm", se = TRUE) + geom_rug(col, alpha = 0.1)
    plot_top <- plot_top + geom_density(col, alpha=.5) 
    plot_right <- plot_right + geom_density(col, alpha=.5)
  }  else {
    scatter <- scatter + geom_point(alpha=.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1)
    plot_top <- plot_top + geom_density(alpha=.5) 
    plot_right <- plot_right + geom_density(alpha=.5)
  }
  
  if (separateBy != "none") {
    scatter <- scatter + theme(legend.position = "none") + facet_wrap(as.formula(paste("~", separateBy)))
  }  else {
    scatter <- scatter
  }
  
  if (separateBy == "none") {
    #arrange the plots together, with appropriate height and width for each row and column
    arrangeGrob(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(3, 1), heights=c(1.5, 3))
  } else {
    scatter
  } 
}

########################################
############## myCorrTest ##############
########################################
# Use to generate summary data for the correlation analysis
myCorrTest <- function (df, gene1, gene2, histo = "All", subtype = "All", colorBy = "none", separateBy = "none",...) {
  if (histo != "All") {
    df <- filter (df, Histology == histo)
  } 
  if (histo == "GBM" & subtype != "All") {
    df <- filter (df, Subtype == subtype)
  }
  if (separateBy == "Histology") {
    cor <- substitute(df %>%
                        group_by(Histology)%>%
                        summarise(n = n(),
                                  r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "Subtype") {
    cor <- substitute(df %>%
                        group_by(Subtype)%>%
                        summarise(n = n(),
                                  r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "none"){
    cor <- cor.test(df[ ,gene1], df[ ,gene2], use = "complete.obs")
#     cor <- data.frame(n = cor$parameter[[1]]+2, r = cor$estimate, p.value = cor$p.value)
  }
  cor
}

########################################
############## pairs panels ############
########################################
myPairsPlot <- function(df,...) { 
  panel.cor <- function(x, y, digits = 3, prefix="", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r)) / 2)
  } 
  panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
  }
  panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = .8, col.smooth = "black", ...) {
    points(x, y, pch = pch, col=rgb(0, 0, 0, 0.5), bg = bg, cex = cex)
    abline(stats::lm(y ~ x), col = "red", ...)
  }
  #   ggpairs(df)
  pairs (df,upper.panel = panel.cor,
         diag.panel = panel.hist,
         lower.panel = panel.lm, pch= 20)
}

######################################
############## kmPlot  ###############
######################################
# Use to plot survival curves, getting the cutoff from the interactive HR plot
kmPlot <- function (cutoff,surv){
  sFit <- survfit(surv)
  sTable <- data.frame(summary(sFit)$table)
  sDiff.log <- survdiff(surv)
  sDiff.mcox <- survdiff(surv,rho = 1)
  log.rank.p <- round(1 - pchisq(sDiff.log$chi, df = 1), 4)
  mantle.cox.p <- round(1 - pchisq(sDiff.mcox$chi, df = 1), 4)
  star.log <- starmaker(log.rank.p)
  star.mcox <- starmaker(mantle.cox.p)
  model <- summary(coxph(surv))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  smax <- max(sFit$time,na.rm=TRUE) 
  smax <- smax-(20*smax)/100
  plot(sFit, yscale = 100 , col = c("red", "blue"), mark.time=FALSE)
  title (main="Kaplan Meier Survival Estimates", xlab="Survival Time (Months)", 
         ylab="% Surviving", font.main = 1, cex.main = 1)
  legend("topright", c(paste("High expr. ", paste(" (n=", sTable$records[1]),", events=",sTable$events[1],", median=",sTable$median[1],")", sep = ""), 
                       paste("Low expr. ", paste(" (n=", sTable$records[2]),", events=",sTable$events[2],", median=",sTable$median[2],")", sep = "")), 
         col= c("red", "blue"), lty = 1, cex = 1)
  text (smax-10, 0.725, paste("HR = ",HR, " (", HR.lower, "-", HR.upper,")", sep=""), cex = 1)
  text (smax-10, 0.65, paste (star.log, "Log-rank p value=", log.rank.p), cex = 1)
  text (smax-10, 0.575, paste (star.mcox, "Wilcoxon p value=", mantle.cox.p), cex = 1)
}


############################################################################
############## Help popup (https://gist.github.com/jcheng5/5913297)  #######
############################################################################
## https://groups.google.com/forum/#!searchin/shiny-discuss/helpPopup/shiny-discuss/ZAkBsL5QwB4/vnmbT47uY7gJ
helpPopup <- function(title, content,
                      placement=c("right", "top", "left", "bottom"),
                      trigger=c("click", "hover", "focus", "manual"),
                      glue = NULL) {
  
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#",
      # class = "btn btn-default", 
      `data-toggle` = "popover",
      title = title, 
      `data-content` = content,
      `data-html` = TRUE, 
      `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],  
      glue, 
      # tags$i(class="icon-info-sign")
      tags$i(class="icon-question-sign")
    ),
    tags$style(type='text/css', ".popover { width: 1200px; relative; top: 20px; left: 20px !important; }")
  )
}

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
maj <- function(InVec) {
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
plotPurity <- function (estimate.df, sample, platform = c("affymetrix","agilent", "illumina")) {
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
