#########################################
##############  ggboxPlot  ##############
#########################################
ggboxPlot <- function(data,  scale = FALSE, xlabel, stat = FALSE, colBox = FALSE, colStrip = FALSE, colorPoints, bw = FALSE, ...) {
  if (scale) {
    ylab <- "Normalized mRNA expression"
  } else {
    ylab <- "mRNA expression (log2)"
  }
  if (colBox) {
    box <- geom_boxplot(aes(fill = group), outlier.size = 0) # It works but not the right way to approach this issue
  } else {
    box <- geom_boxplot(outlier.size = 0)
  }
  if (colStrip) {
    col <- aes_string(color = colorPoints)
    strip <- geom_jitter(position = position_jitter(width = .2), col, size = 2, alpha = 0.75)
  } else {
    strip <- geom_jitter(position = position_jitter(width = .2), size = 2, alpha = 0.5)
  }
  p <- ggplot(data, aes(x=group, y = mRNA)) + ylab(ylab) + xlab(xlabel) 
  p <- p + box + strip 
  if (bw) {
    p <- p + theme_bw () 
  }
  if (stat) {
    tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
    tukey <<- tukey ##  see scoping rules http://shiny.rstudio.com/articles/scoping.html
    tukey$Significance <- as.factor(starmaker(tukey$p.adj,p.levels=c(.001, .01, .05, 1), 
                                              symbols=c("***", "**", "*", "ns")))
    tukey$comparison <- row.names(tukey)
    
    t <- ggplot(tukey, aes(reorder(comparison, diff), diff, ymin = lwr, ymax= upr, colour = Significance)) +
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
}

######################################
############## Get HR  ###############
######################################
getHR <- function (df, gene) {
  mRNA <- df[ ,gene]
  surv.status <- df[ ,"status"]
  surv.time <- df[ ,"survival"]
  my.Surv <- Surv(surv.time, surv.status == 1)
  mRNA.values <- mRNA[!is.na(mRNA)]
  # Generate a vector of continuos values, excluding the first an last value
  mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) 
  scan.surv <-function(i, conf.level=95) {
    log.rank <- survdiff(my.Surv ~ mRNA <= i, data = df, rho = 0)
    model <- summary(coxph(my.Surv ~ mRNA <= i))
    HR <- model$conf.int[1]
    HR.lower <- model$conf.int[3]
    HR.upper <- model$conf.int[4]
    data.frame(i, HR, HR.lower, HR.upper, log.rank$obs [1], log.rank$obs [2])
  }
  HRdata <- data.frame (t(sapply(mRNA.values, scan.surv)))
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
survivalPlot <- function (df, gene, group, cutoff, numeric, subtype, gcimp = FALSE, primary = FALSE) {
  df <- subset(df, !is.na(df$status))
  # Select the samples, for the specified histology,that has survival data
  if (group != "All") {
    df <- subset (df, Histology == group) 
  }
  # For GBM, select only primary tumors
  if (primary & any(!is.na(df$Recurrence))) {
    df <- subset (df, Recurrence == "Primary")
  }
  # Select a specific subtype
  if (group == "GBM" & subtype != "All") {
    df <- subset (df, Subtype == subtype)
  }
  # Remove G-CIMP, when selected
  if (gcimp){
    df <- subset (df, Subtype != "G-CIMP")
  } 
  df <- subset(df, !is.na(df$status))
  mRNA <- df[ ,gene]
  surv.status <- df[ ,"status"]
  surv.time <- df[ ,"survival"]
  my.Surv <- Surv(time = surv.time, event = surv.status== 1)
  smax <- max(surv.time, na.rm = TRUE)
  tmax <- smax-(25*smax)/100
  mRNA.q <- quantile(mRNA, probs=c(0.25, 0.5, 0.75), na.rm = TRUE)
  
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
    cut <- switch(cutoff, 
           "median" = mRNA.q[2],
           "lower quartile" = mRNA.q [1],
           "upper quartile" = mRNA.q [3],
           "Use a specific mRNA value" = numeric)

    f <- function(x) ifelse(x >= cut, c("high"),c("low"))
    strat <- f(mRNA)
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
         main = main,
         col = c("red", "blue"), 
         mark.time = FALSE)
    legend("topright", legend = c(paste(gene," High ", paste("(n=", surv$records[1]),", events=", surv$events[1],", median=",surv$median[1],")", sep = ""), 
                                  paste(gene," Low ", paste("(n=", surv$records[2]),", events=", surv$events[2],", median=",surv$median[2],")", sep = "")),
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
           legend = c(paste("1st ",paste("(n=", z$records[1]),", median=",z$median[1],")", sep = ""),
                      paste("2nd ",paste("(n=", z$records[2]),", median=",z$median[2],")", sep = ""),
                      paste("3rd ",paste("(n=", z$records[3]),", median=",z$median[3],")", sep = ""),
                      paste("4th ",paste("(n=", z$records[4]),", median=",z$median[4],")", sep = "")), 
           col= c(1:4), lty=1, cex=1)
  }
  
}

#####################
## Get correlations ##
#####################
# To use to geet correlation data (r an p value) on the fly. 
# using Hadley suggestion: https://stat.ethz.ch/pipermail/r-help/2008-November/181049.html
getCorr <- function (df, gene, histology, corrMethod) {
  if (histology != "All") {
    df <- subset (df, Histology == histology)
  } else {
    df <- df
  }
  df <- df[,8:ncol(df)]
  mRNA <- df[ ,gene, drop = F]
  r <- apply(mRNA, 2, function(x) { apply(df, 2, function(y) { cor(x,y, method = corrMethod) })})
  df <- nrow(mRNA) - 2
  t <- sqrt(df) * r / sqrt(1 - r ^ 2)
  p <- pt(t, df)
  p <- 2 * pmin(p, 1 - p)
  corr <- data.frame(row.names(r), round(r,3), round(p,10))
  names(corr) <- c("Gene","r","p.value")
  corr
}

######################################################
############## 2 genes correlation plot ##############
######################################################
myCorggPlot <- function (df, gene1, gene2, histo = "All", subtype = "All", colorBy = "none", separateBy = "none",...) {
  if (histo != "All") {
    df <- subset (df, Histology == histo)
  } else {
    df <- df
  }
  if (histo == "GBM" & subtype != "All") {
    df <- subset (df, Subtype == subtype)
  }
  #  empy plot to used in grid.arrange 
  empty <- ggplot() + geom_point(aes(1,1), colour="white") + 
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
  # marginal density of x - plot on top
  plot_top <- ggplot(df, mapping = aes_top) + theme(legend.position = "none", axis.title.x = element_blank())
  # marginal density of y - plot on the right
  plot_right <- ggplot(df, mapping = aes_right) + coord_flip() + theme(legend.position = "none",axis.title.y = element_blank())
  
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
    cor <- substitute(df %>%
                        group_by(Histology)%>%
                        summarise(r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "Subtype") {
    cor <- substitute(df %>%
                        group_by(Subtype)%>%
                        summarise(r = cor.test(x, y, use = "complete.obs")$estimate,
                                  p.value = cor.test(x, y, use = "complete.obs")$p.value), 
                      list(x = as.name(gene1), y = as.name(gene2)))
    cor <- data.frame(eval(cor))
  } else if (separateBy == "none"){
    cor <- cor.test(Gene1, Gene2, use = "complete.obs")
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
  
  
  # splits <- str_c(as.character(mapping$group), as.character(mapping$colour), sep = ", ", collapse = ", ")
  # splits <- str_c(colorCol, sep = ", ", collapse = ", ")
  final_text <- ""
  if(length(colorCol) < 1)
    colorCol <- "ggally_NO_EXIST"
  # browser()
  if(colorCol != "ggally_NO_EXIST" && colorCol %in% colnames(data)) {
    
    txt <- str_c("ddply(data, .(", colorCol, "), summarize, ggally_cor = cor(", xCol,", ", yCol,"))[,c('", colorCol, "', 'ggally_cor')]")
    
    con <- textConnection(txt)
    on.exit(close(con))
    cord <- eval(parse(con))
    
    # browser()
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
    # print(order(ord[ord >= 0]))
    # print(lev)
    cord <- cord[order(ord[ord >= 0]), ]
    
    cord$label <- str_c(cord[[colorCol]], ": ", cord$ggally_cor)
    
    # calculate variable ranges so the gridlines line up
    xmin <- min(xVal)
    xmax <- max(xVal)
    xrange <- c(xmin-.01*(xmax-xmin),xmax+.01*(xmax-xmin))
    ymin <- min(yVal)
    ymax <- max(yVal)
    yrange <- c(ymin-.01*(ymax-ymin),ymax+.01*(ymax-ymin))
    
    
    # print(cord)
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
    # print(range(yVal))
    # print(yPos)
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