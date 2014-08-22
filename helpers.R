############################################
############## my ggboxPlot  ##############
############################################
ggboxPlot <- function(exprs, cna, gene, plotType, scale = FALSE, stat = FALSE, colBox = FALSE, colStrip = FALSE, bw = FALSE, ...) {
  if (!gene%in%names(exprs)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
   mRNA <- exprs[ ,gene] 
  if (scale) {
    mRNA <- scale(mRNA)
    ylab <- "Normalized mRNA expression"
  } else if (!scale){
    mRNA <- mRNA
    ylab <- "mRNA expression (log2)"
  }
  if (plotType == "Copy number") {
    group <- cna[, gene]
    group <- factor(group, levels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp")) # I should have done this directly on all the cna in the Rds file
    group <- droplevels(group)
  } else {
    group <- exprs[ ,plotType]
  }
  if (any(!is.na(group))) {
  data <- data.frame(mRNA,group)
  data <- na.omit(data)
  n.class <- length(levels (group)) #  Not used yet
  if (colBox) {
    box <- geom_boxplot(aes(fill = group), outlier.size = 0) # It works but not the right way to approach this issue
  } else {
    box <- geom_boxplot(outlier.size = 0)
  }
  if (colStrip) {
    strip <- geom_jitter(position = position_jitter(width = .2), 
                         aes(colour = group), size = 2, alpha = 0.75)
#       theme(legend.position = "none")
  } else {
    strip <- geom_jitter(position = position_jitter(width = .2), size = 2, alpha = 0.5)
  }
  p <- ggplot(data, aes(x=group, y = mRNA)) + ylab(ylab) + xlab("")
  p <- p + box + strip
  if (bw) {
    p <- p + theme_bw ()
  }
  if (stat) {
    tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
    tukey <<- tukey ##  see scoping rules http://shiny.rstudio.com/articles/scoping.html
    tukey$Significance <- as.factor(starmaker(tukey$p.adj,p.levels=c(.001, .01, .05, 1), 
                                              symbols=c("***", "**", "*", "ns")))
    
    t <- ggplot(tukey, aes(row.names(tukey), diff, ymin = lwr, ymax= upr, colour = Significance)) +
      geom_point() + geom_errorbar(width = 0.25) + 
      ylab("Differences in mean levels") + xlab("") + 
      geom_hline(xintercept = 0, colour="darkgray", linetype = "longdash") + coord_flip()
    if (bw) {
      t <- t + theme_bw ()
    }
    grid.arrange(p, t, ncol=2, widths = c(3,2))
  } else {
    print(p) 
  }
  } else {
    stop (paste(plotType, "not available for this dataset"))
  }
}

############
## Get HR ##
############
getHR <- function (df, gene, gcimp = FALSE) {
  df <- subset (df, Histology == "GBM")
  if (gcimp){
    df <- subset (df, Subtype != "G-CIMP")
  }
  if (!gene%in%names(df)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
  mRNA <- df[ ,gene]
  surv.status <- df[ ,"status"]
  surv.time <- df[ ,"survival"]
  my.Surv <- Surv(surv.time, surv.status == 1)
  mRNA.values <- mRNA[which(surv.status != "NA")] # mRNA expression in patient with survival data
  mRNA.values <- mRNA.values[!is.na(mRNA.values)]
  mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) # Generate a a vector of continuos values, excluding the first an last value
  scan.surv <-function(i,conf.level=95) {
    log.rank <- survdiff(my.Surv ~ mRNA <= i, data = df, rho = 0)
    model <- summary(coxph(my.Surv ~ mRNA <= i))
    HR <- model$conf.int[1]
    HR.lower <- model$conf.int[3]
    HR.upper <- model$conf.int[4]
    data.frame(i, HR, HR.lower, HR.upper, log.rank$obs [1], log.rank$obs [2])
  }
  HRdata <- data.frame (t(sapply(mRNA.values, scan.surv)))
  HRdata <- data.frame (sapply(HRdata,unlist))
  HRdata <- subset(HRdata, HRdata[,5] >= 15 & HRdata[,6] >= 15) # Exclude groups with less than 15 samples. ARBITRARY
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
  points(HRdata[, "i"], rep(ymin - 0.4, length(HRdata[, "i"])), pch="|")
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
survivalPlot <- function (df, gene, group, cutoff, subtype, gcimp = FALSE) {
  if (!gene%in%names(df)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
  df <- subset (df, Histology == group & Histology != "Non-tumor") 
  if (group == "Non-tumor") {
    stop ("Sorry no survival data are available for this group") ## Need to throw an error aslo if  SEE REMBRANDT PROBLEM
  }
  if (group == "GBM" & any(!is.na(df$Recurrence))) {
    df <- subset (df, Histology == "GBM" & Recurrence == "Primary")
  } 
  if (subtype != "All") {
    df <- subset (df, Subtype == subtype)
  }
  if (gcimp){
    df <- subset (df, Subtype != "G-CIMP")
  }
  mRNA <- df[ ,gene]
  surv.status <-df[ ,"status"]
  surv.time <- df[ ,"survival"]
  smax <- max(surv.time, na.rm=TRUE)
  tmax <- smax-(25*smax)/100
  my.Surv <- Surv(surv.time, surv.status== 1)
  mRNA.q <- quantile(mRNA, probs=c(0.25, 0.5, 0.75), na.rm = TRUE)
  main <- paste("Histology: ", group, 
                "; Subtype: ", subtype,
                "; Cutoff: ", cutoff, sep = "")
  if (cutoff != "quartiles") {
    if (cutoff == "median") {
      cut <- mRNA.q[2]
    } else if (cutoff == "lower quartile") {
      cut <- mRNA.q [1]
    } else if (cutoff == "upper quartile") {
      cut <- mRNA.q [3]
    }
    f<-function(x) ifelse(x >= cut, c("high"),c("low"))
    mRNA <- f(mRNA)
    expr.surv <- survfit(my.Surv ~ mRNA, data = df, conf.type = "none")
    log.rank <- survdiff(my.Surv ~ mRNA, data = df, rho = 0)
    mantle.cox <- survdiff(my.Surv~ mRNA, data = df, rho = 1)
    surv <- data.frame(summary(expr.surv)$table)
    model <- summary(coxph(my.Surv ~ mRNA, data = df))
    HR <- round(model$conf.int[1],2)
    HR.lower <- round(model$conf.int[3],2)
    HR.upper <- round(model$conf.int[4],2)
    log.rank.p <- round(1 - pchisq(log.rank$chi, df = 1), 4)
    mantle.cox.p <- round(1 - pchisq(mantle.cox$chi, df = 1), 4)
    star.log <- starmaker(log.rank.p)
    star.mcox <- starmaker(mantle.cox.p)
    plot(expr.surv, xlab = "Survival time (Months)", ylab = "% Surviving", yscale = 100, xlim = c(0,smax),
         main = main,
         col = c("red", "blue"), 
         mark.time = FALSE)
    legend("topright", legend = c(paste(gene," High ", paste("(n=", surv$events[1]),", median=",surv$median[1],")", sep = ""), 
                                  paste(gene," Low ", paste("(n=", surv$events[2]),", median=",surv$median[2],")", sep = "")),
           col= c("red", "blue"), lty = 1, cex = 1)
    text (tmax-10, 0.725, paste("HR = ",HR, " (", HR.lower, "-", HR.upper,")", sep=""), cex = 1)
    text (tmax-10, 0.65, paste (star.log, "Log-rank p value=", log.rank.p), cex = 1)
    text (tmax-10, 0.575, paste (star.mcox, "Wilcoxon p value=", mantle.cox.p), cex = 1)
  }
  
  if (cutoff == "quartiles"){
    exprcat <- cut(mRNA, quantile(mRNA,na.rm = T), data=df)
    expr.surv <- survfit(my.Surv ~ strata(exprcat), data=df, conf.type="none")
    z <- data.frame(summary(expr.surv)$table) 
    plot(expr.surv, xlab="Months", ylab="% Surviving", yscale = 100, xlim = c(0,smax), 
         main = main,
         col= c(1:4),
         mark.time=FALSE)
    legend("topright", title = "Quantile", 
           legend = c(paste("1st ",paste("(n=", z$events[1]),", median=",z$median[1],")", sep = ""),
                      paste("2nd ",paste("(n=", z$events[2]),", median=",z$median[2],")", sep = ""),
                      paste("3rd ",paste("(n=", z$events[3]),", median=",z$median[3],")", sep = ""),
                      paste("4th ",paste("(n=", z$events[4]),", median=",z$median[4],")", sep = "")), 
           col= c(1:4), lty=1, cex=1)
  }
  
}

#####################
## Get correlations ##
#####################
# To use to geet correlation data (r an p value) on the fly. TOO SLOW
getCorr <- function (df, gene, histology) {
  if (!gene%in%names(df)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
  if (histology != "All") {
    df <- subset (df, Histology == histology)
  } else {
    df <- df
  }
  df <- df[,8:ncol(df)]
  mRNA <- df[ ,gene]
  corr <- NULL
  for (Gene in names(df)){
    y <- df[ ,Gene]
    cor <- cor.test (mRNA, y, use="complete.obs")
    r <- round(cor$estimate, digits =  3)
    p <- round(cor$p.value, digits = 10)
    corr <- rbind(corr, data.frame(Gene, r, p))
  }
  corr
}

######################################################
############## 2 genes correlation plot ##############
######################################################
myCorggPlot <- function (df, gene1, gene2, histo = "All", subtype = "All", colorBy = "none", separateBy = "none",...) {
  if (!gene1%in%names(df) | !gene2%in%names(df)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
  if (histo != "All") {
    df <- subset (df, Histology == histo)
  } else {
    df <- df
  }
  if (histo == "GBM" & subtype != "All") {
    df <- subset (df, Subtype == subtype)
  }
  #  empy plot to used in grid.arrange 
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(plot.background = element_blank(), panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), panel.border = element_blank(), 
      panel.background = element_blank(), axis.title.x = element_blank(),
      axis.title.y = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank()
    )
  
  aes_scatter <- aes_string(x = gene1, y = gene2)
  aes_top <- aes_string(x = gene1)
  aes_right <- aes_string(x = gene2)
  # scatterplot of x and y variables
  scatter <- ggplot(df,mapping = aes_scatter) + theme(legend.position=c(1,1),legend.justification=c(1,1))
  #marginal density of x - plot on top
  plot_top <- ggplot(df, mapping = aes_top) + theme(legend.position = "none",axis.title.x = element_blank())
  #marginal density of y - plot on the right
  plot_right <- ggplot(df, mapping = aes_right) + coord_flip() + theme(legend.position = "none",axis.title.y = element_blank())
  
  if (colorBy == "Histology") {
    scatter <- scatter + 
      geom_point(aes(color = Histology), alpha=.5) + 
      geom_smooth(aes(color = Histology), method = "lm", se = TRUE)
    plot_top <- plot_top + geom_density(aes(color = Histology), alpha=.5) 
    plot_right <- plot_right + geom_density(aes(color = Histology), alpha=.5)
  } else if (colorBy == "Subtype") {
    scatter <- scatter + 
      geom_point(aes(color = Subtype),alpha=.5) + 
      geom_smooth(aes(color = Subtype), method = "lm", se = TRUE)
    plot_top <- plot_top + geom_density(aes(color = Subtype), alpha=.5) 
    plot_right <- plot_right + geom_density(aes(color = Subtype), alpha=.5)
  } else {
    scatter <- scatter + 
      geom_point(alpha=.5) + 
      geom_smooth(method = "lm", se = TRUE)
    plot_top <- plot_top + geom_density(alpha=.5) 
    plot_right <- plot_right + geom_density(alpha=.5)
  }
  
  if (separateBy == "Histology") {
    scatter <- ggplot(df,mapping = aes_scatter) + 
      geom_point(aes(color = Histology), alpha=.5) + 
      geom_smooth(aes(color = Histology), method = "lm", se = TRUE) +
      facet_wrap (~ Histology)
  } else if (separateBy == "Subtype") {
    scatter <- ggplot(df,mapping = aes_scatter) + 
      geom_point(aes(color = Subtype),alpha=.5) + 
      geom_smooth(aes(color = Subtype), method = "lm", se = TRUE) +
      facet_wrap (~ Subtype)  
  } else {
    scatter <- scatter
  }
  
  if (separateBy == "none") {
  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(3, 1), heights=c(1.5, 3))
  } else {
    print(scatter) 
  }

}

########################################
############## myCorrTest ##############
########################################
# Use to generate summary data for the correlation analysis
myCorrTest <- function (df, gene1, gene2, histo = "All", subtype = "All", colorBy = "none", separateBy = "none",...) {
  if (!gene1%in%names(df) | !gene2%in%names(df)) {
    stop ("Incorrect gene entry or gene not available for this dataset")
  }
  if (histo != "All") {
    df <- subset (df, Histology == histo)
  } else {
    df <- df
  }
  if (histo == "GBM" & subtype != "All") {
    df <- subset (df, Subtype == subtype)
  }
  Gene1 <- df[ ,gene1]
  Gene2 <- df[ ,gene2]
  if (separateBy == "Histology") {
    cor <- substitute(df %.%
                        group_by(Histology)%.%
                        summarise(r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "Subtype") {
    cor <- substitute(df %.%
                        group_by(Subtype)%.%
                        summarise(r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "none"){
  cor <- cor.test(Gene1, Gene2, use = "complete.obs")
  }
  cor
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
  model <- summary(coxph(surv))
  HR <- round(model$conf.int[1],2)
  HR.lower <- round(model$conf.int[3],2)
  HR.upper <- round(model$conf.int[4],2)
  smax <- max(sFit$time,na.rm=TRUE) 
  smax <- smax-(20*smax)/100
  plot(sFit, yscale = 100 , col = c("red", "blue"), mark.time=FALSE)
  title (main="Kaplan Meier Survival Estimates", xlab="Survival Time (Months)", 
         ylab="% Surviving", font.main = 1, cex.main = 1)
  legend("topright", c(paste("High expr. ", paste(" (n=", sTable$events[1]),", median=",sTable$median[1],")", sep = ""), 
                       paste("Low expr. ", paste(" (n=", sTable$events[2]),", median=",sTable$median[2],")", sep = "")), 
         col= c("red", "blue"), lty = 1, cex = 1)
  text (smax-10, 0.725, paste("HR = ",HR, " (", HR.lower, "-", HR.upper,")", sep=""), cex = 1)
  text (smax-10, 0.65, paste ("Log-rank p value=", log.rank.p), cex = 1)
  text (smax-10, 0.575, paste ("Wilcoxon p value=", mantle.cox.p), cex = 1)
}


######################################
############## getData  #############
######################################
getData <- function (df, gene){
  data <- df[ ,c(names(df)[1:7], gene)]
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

############################################################################
############## binding to a bootstrap modal, https://github.com/mostly-harmless/radiant/commit/9b554532e417b4a8c0e8a43f16111f9b57c343d4#diff-848f402b1e98cb59e92553dac731963dR161  ###########
############################################################################
helpModal <- function(title, link, content) {
  html <- sprintf("<div id='%s' class='modal hide fade in' style='display: none; '>
                  <div class='modal-header'><a class='close' data-dismiss='modal' href='#'>&times;</a>
                  <h3>%s</h3>
                  </div>
                  <div class='modal-body'>%s</div>
                  </div>
                  <a title='Help' data-toggle='modal' href='#%s' class='icon-question-sign'></a>", link, title, content, link)
  Encoding(html) <- 'UTF-8'
  HTML(html)
}

