library(shiny)
library(survival)
library(weights)
library(ggplot2)
library(gridExtra)
library(rCharts)
# library(shinysky)
library(dplyr)
library(GSVA)

source("global.R")
`%then%` <- shiny:::`%OR%`

#######################################
############## Datasets  ##############
#######################################
gene_names <- readRDS("data/gene_names.Rds")
plotList <- list("TCGA GBM" = c("Histology", "Copy number", "Subtype","Recurrence"),
                 "TCGA Lgg" = c("Histology", "Grade", "Copy number", "Subtype"),
                 "Rembrandt" = c("Histology", "Grade", "Subtype", "Recurrence"),
                 "Gravendeel" = c("Histology", "Grade", "Subtype"),
                 "Phillips" = c("Grade", "Subtype", "Recurrence"),
                 "Murat" = c("Histology", "Subtype", "Recurrence"),
                 "Freije" = c("Histology", "Grade", "Subtype"),
                 "Reifenberger" = c("Subtype"),
                 "Bao" = c("Histology", "Subtype", "Recurrence"),
                 "Gill" = c("Histology", "Subtype"),
                 "Gorovets" = c("Histology", "Grade", "Subtype"),
                 "Nutt" = c("Histology", "Subtype"))
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
subtype_list <- readRDS("data/subtype_list.Rds")

#######################################
############## server.R  ##############
#######################################
shinyServer(
  function(input, output, session) {
    
    options(shiny.maxRequestSize=30*1024^2)    
    
    #' Return the requested dataset
    datasetInput <- reactive({
      switch(input$dataset, 
             "TCGA GBM" = gbm.tcga,
             "TCGA Lgg" = lgg.tcga,
             "Rembrandt" = rembrandt,
             "Gravendeel" = gravendeel,
             "Phillips" = phillips,
             "Murat" = murat,
             "Freije" = freije,
             "Reifenberger" = reifenberger,
             "Bao" = bao,
             "Gill" = gill,
             "Gorovets" = gorovets,
             "Nutt" = nutt)
    })
    
    #' Switch the datset for the correlation
    datasetInputCor <- reactive({
      switch(input$datasetCor, 
             "TCGA GBM" = gbm.tcga,
             "TCGA Lgg" = lgg.tcga,
             "Rembrandt" = rembrandt,
             "Gravendeel" = gravendeel,
             "Phillips" = phillips,
             "Murat" = murat,
             "Freije" = freije,
             "Reifenberger" = reifenberger,
             "Bao" = bao,
             "Gill" = gill,
             "Gorovets"= gorovets,
             "Nutt" = nutt)
    })
    
    #' Expression data
    exprs <- reactive({
      datasetInput()[["expr"]]
    })
    
    #' Phenotype data
    pDatas <- reactive({
      datasetInput()[["pData"]]
    })
    
    #' CNA data
    cnas <- reactive({
      datasetInput()[["cna"]]
    })
    
    #' Text matching with the gene names list
    updateSelectizeInput(session, inputId = "gene", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "geneCor", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "gene2", choices = gene_names, server = TRUE) 
    
    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plotSelected <- reactive ({
      if (input$plotTypeSel %in% plotList[[input$dataset]]){ 
        input$plotTypeSel
      } else {
        NULL
      }
    })
    
    #' Change the plot type available for a specific dataset
    observe({
      updateSelectInput(session, inputId = "plotTypeSel", choices = plotList[[input$dataset]], selected = plotSelected()) 
    }, priority = 10)
    
    # Caption with gene and dataset
    output$caption <- renderText({
      if (input$gene == "" )
        return()
      title <- paste(input$gene, "in", input$dataset, "dataset")
    })
    
    #' Return the available histology, to be used in the updateSelectInput for survival correlation 
    histo <- reactive({
      levels(datasetInput()[["pData"]][,"Histology"])
    })
    
    #' When switching datasets for surv, if the selected histo is not available it will choose GBM (the last histo of the list)
    histoSurvSelected <- reactive ({
      if (input$histologySurv %in% histo()){
        input$histologySurv
      } else {
        tail(histo(), n=1)
      }
    })
    
    #' When switching datasets for correlation plots, if the selected histo is not available it will choose "All"
    histoCorrSelected <- reactive ({
      if (input$histologyCorr %in% histo()){
        input$histologyCorr
      } else {
        return("All")
      }
    })
    
    observe({
      # This will change the value of input$histologySurv, based on histological group available for that dataset  
      updateSelectInput(session, inputId = "histologySurv", choices = histo(), selected = histoSurvSelected())
      # This will change the value of input$histologyCorr, based on histological group available for that dataset
      updateSelectInput(session, inputId = "histologyCorr", choices = c("All", histo()), selected = histoCorrSelected())
    })
      
    #' Return the available histology, to be used in the updateSelectInput for correlation table
    histoCor <- reactive({
      levels(datasetInputCor()[["pData"]][,"Histology"])
    })
    
    #' Change the value of input$histologyCorrTable, based on histological group available for that dataset
    observe({
      updateSelectInput(session, inputId = "histologyCorrTable", choices = c("All", histoCor()), selected = "All")
    })
        
    #' Create the dataframe to call in ggbox,Tukey and ttest
    data <- reactive({     
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the plotType is not available.
        need(plotSelected() %in% plotList[[input$dataset]],"")
      ) 
      mRNA <- exprs()[ ,input$gene]
      if (input$scale) {
        mRNA <- scale(mRNA)
      }
      if (plotSelected() == "Copy number") {
        validate(need(input$gene %in% names(cnas()), "Copy number not available for this gene in this dataset"))
        group <- cnas()[ ,input$gene]
        group <- factor(group, levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        group <- droplevels(group)
      } else {
        group <- exprs()[ ,plotSelected()]
      }
      data <- data.frame(mRNA, group, pDatas()[,-c(1,6,7)]) # To exclude sample name and survival data
      data <- subset(data, !is.na(group))
      data
    })
    
    #' Generate radiobuttons for the various categories in the pData
    output$colorPoints <- renderUI({
      if (input$gene == "")
        return()
      data <- data()[,-c(1:2)] # to remove mRNA and group
      data <- data[,colSums(is.na(data)) < nrow(data)] # to remove NA columns 
      colnames <- names(data)
      # Create the radiobuttons for the different pData categories
      radioButtons("colorP", "Color by:", choices  = colnames, selected = plotSelected())
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({
      # Trying to avoid an error when switching datasets in case the colStrip is not available.
      if(input$colStrip){
        validate(need(input$colorP %in% names(data()[,-c(1:2)]),""))
      }
      # I needed to create the ggboxPlot function (see helper file) to use it in the output$downloadPlot
      ggboxPlot(data = data (), scale = input$scale, stat = input$stat, colBox = input$colBox, 
                colStrip = input$colStrip, colorPoints = input$colorP, bw = input$bw)  
    })
    
    #' Tukey post-hoc test
    output$tukeyTest <- renderPrint(width = 800, {    
      data <- data()
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff), ]
      tukey
    })
    
    #' Pairwise t test
    output$pairwiseTtest <- renderPrint({     
      data <- data()
      pttest <- pairwise.t.test(data$mRNA, data$group, na.rm= TRUE, p.adj = "bonferroni", paired = FALSE)[[3]]
      pttest
    })
    
    #' Get the selected download file type.
    downloadPlotFileType <- reactive({
      input$downloadPlotFileType  
    })    
    observe({
      plotFileType    <- input$downloadPlotFileType
      plotFileTypePDF <- plotFileType == "pdf"
      plotUnit    <- ifelse(plotFileTypePDF, "inches", "pixels")
      plotUnitDef <- ifelse(plotFileTypePDF, 7, 480)
      
      updateNumericInput(
        session,
        inputId = "downloadPlotHeight",
        label = sprintf("Height (%s)", plotUnit),
        value = plotUnitDef)
      
      updateNumericInput(
        session,
        inputId = "downloadPlotWidth",
        label = sprintf("Width (%s)", plotUnit),
        value = plotUnitDef)
    })
    
    #' Get the download dimensions.
    downloadPlotHeight <- reactive({
      input$downloadPlotHeight
    })
    downloadPlotWidth <- reactive({
      input$downloadPlotWidth
    })
    
    #' Download the Plot
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_", input$plotType, ".", downloadPlotFileType(), sep = "")
      },     
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(con) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(con, width = downloadPlotWidth(), height = downloadPlotHeight())
        ggboxPlot(data = data (), scale = input$scale, stat = input$stat, colBox = input$colBox, 
                  colStrip = input$colStrip, colorPoints = input$colorP, bw = input$bw)  
        dev.off(which=dev.cur())
      }
    )
    
    # Need a wrapper around the hrClick input so we can manage whether or 
    # not the click occured on the current Gene. If it occured on a previous
    # gene, we'll want to mark that click as 'stale' so we don't try to use it
    # later. https://gist.github.com/trestletech/5929598
    currentClick <- list(click = NULL, stale = FALSE)
    
    handleClick <- observe({
      if (!is.null(input$hrClick) && !is.null(input$hrClick$x)){
        currentClick$click <<- input$hrClick
        currentClick$stale <<- FALSE
      }
    }, priority=100)
    
    #' Extract the GBM samples for the interactive HR plot.
    survData <- reactive({    
      df <- exprs()
      df <- subset(df, !is.na(df$status))
      # subset to primary GBM, in case both primary and recurrent samples are available
      if (any(!is.na(df$Recurrence))) {
        df <- subset (df, Histology == "GBM" & Recurrence == "Primary")
      }
      df <- subset (df, Histology == "GBM")
      # exclude G-CIMP is selected
      if (input$gcimpSurv){
        df <- subset (df, Subtype != "G-CIMP")
      }
      df
    })
    
    #' Extract the GBM expression values for the interactive HR plot.
    geneExp <- reactive({
      geneExp <- survData()[ ,input$gene]
      currentClick$stale <<- TRUE
      geneExp
    })
    
    #' Generate the cutoff value for the interactive HR plot.
    getCutoff <- reactive({
      input$hrClick
      geneExp()        
      # See if there's been a click since the last gene change.
      if (!is.null(currentClick$click) && !currentClick$stale){
        return(currentClick$click$x)
      }       
      median(geneExp())
    })
    
    #' Extract the Hazard ratios for the input gene.
    HR <- reactive ({
      HR <- getHR(survData(), input$gene)
    })
    
    #' Requirements for all the survival plots
    survNeed <- reactive({
      validate(
        need(input$dataset!= "Bao" & input$dataset!= "Reifenberger" & input$dataset!= "Gill" , "Sorry, no survival data are available for this dataset")%then%
                 need(input$histologySurv != "Non-tumor","Sorry, no survival data are available for this group")%then%
                 need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
                 need(input$gene %in% names(exprs()),"Gene not available for this dataset")
        )
    })
    
    #' Render a plot to show the the Hazard ratio for the gene's expression values
    output$hazardPlot <- renderPlot({        
      validate(need(histoSurvSelected() == "GBM", "Interactive HR plot currently available only for GBM samples"))
      survNeed()
      # Plot the hazardplot 
      hazardPlot(HR(), input$quantile)
      # Add a vertical line to show where the current cutoff is.
      abline(v = getCutoff(), col = 4)
    }, bg = "transparent")
        
    #' A reactive survival formula
    survivalFml <- reactive({
      # Create the groups based on which samples are above/below the cutoff
      expressionGrp <- as.integer(geneExp() < getCutoff())
      # Create the survival object
      surv <- with(survData(), Surv(survival, status == 1))
      return(surv ~ expressionGrp)
    })
    
    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      validate(need(histoSurvSelected() == "GBM", "Interactive HR plot currently available only for GBM samples"))
      survNeed()
      cutoff <- getCutoff()
      surv <- survivalFml()
      kmPlot(cutoff, surv)
    })
    
    #' Create a slider for the manual cutoff of the Kaplan Meier plot
    mRNAsurv <- reactive({     
      validate(
        need(input$dataset!= "Bao" & input$dataset!= "Reifenberger" & input$dataset!= "Gill", "")%then%
          need(input$gene != "", "")%then%
          need(input$gene %in% names(exprs()),""),
        need(input$histologySurv %in% histo(),""),
        need(input$histologySurv != "Non-tumor","")
      )     
      df <- subset(exprs(), Histology == input$histologySurv)
      if (input$histologySurv == "GBM" & input$subtypeSurv != "All") {
        df <- subset (df, Subtype == input$subtypeSurv)
      }
      if (input$gcimpSurv){
        df <- subset (df, Subtype != "G-CIMP")
      }
      df <- subset(df, !is.na(df$status))
      mRNA <- df[ ,input$gene]
      mRNA.values <- round(mRNA[!is.na(mRNA)],2)
      # Generate a vector of continuos values, excluding the first an last value
      mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) 
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxmRNA <- renderPlot({    
      validate(need(input$mInput, ""))      
      mRNA <- mRNAsurv()
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|")
      # Add a red line to show which  is the current cutoff.
      points(x = input$mInput, y = 0, pch = "|", col="red",cex = 1.5)
    }, bg = "transparent", width = 225)
    
    #' Generate the slider for the manual cutoff
    output$numericCutoff <- renderUI({
      sliderInput(inputId = "mInput",label = NULL, min = min(mRNAsurv()), max = max(mRNAsurv()), value = median(mRNAsurv()))
    })
    
    #' Create a Kaplan Meier plot with cutoff based on quantiles or manual selection
    output$survPlot <- renderPlot({     
      survNeed ()
      validate(need(input$histologySurv %in% histo(),""))   
      # Use try because I need to suppress a message throwed the first time manual cutoff is selected
      try(survivalPlot (exprs(), input$gene, group = input$histologySurv, cutoff = input$cutoff, numeric = input$mInput,
                        subtype = input$subtypeSurv, gcimp = input$gcimpSurv), silent = TRUE) 
    })
    
    #' Download the survPlot
    output$downloadsurvPlot <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_survPlot.pdf", sep = "")
      },      
      content = function(file) {
        pdf(file)
        survivalPlot(exprs(), input$gene, group = input$histologySurv, cutoff = input$cutoff, 
                     subtype = input$subtypeSurv, gcimp = input$gcimpSurv)
        dev.off()
      }
    )
    
    #' Download the kmPlot
    output$downloadkmPlot <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_kmPlot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file)
        kmPlot(getCutoff(), survivalFml())
        dev.off()
      }
    )
    
    #' Generate reactive Inputs for the corrPlot to be used also to download the complete plot
    colorByInput <- reactive({
      switch(input$colorBy, 
             none = "none",
             Histology = "Histology",
             Subtype = "Subtype")
    })
    
    separateByInput <- reactive({
      switch(input$separateBy, 
             none = "none",
             Histology = "Histology",
             Subtype = "Subtype")
    })
    
    #' Generate the correlation plot
    output$corrPlot <-renderPlot({    
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        need(input$gene2 != "", "Please enter Gene 2")%then%
          need(input$gene2 %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the choosen histology is not available.
        need(input$histologySurv %in% histo(),"")
      )
      if (input$dataset == "TCGA Lgg") {
        validate(need(input$colorBy != "Subtype" & input$separateBy != "Subtype", "Subtype available for GBM samples only")) 
      }
      myCorggPlot(exprs(), input$gene, input$gene2, input$histologyCorr, input$subtype, 
                  colorBy = colorByInput(), separateBy = separateByInput())
    })
    
    #' Generate a summary of the correlation test
    output$corrTest <-renderPrint({     
      if (input$gene == "" | input$gene2 == "")
        return()
      validate(
        # Trying to avoid an error when switching datasets in case the choosen histology is not available.
        need(input$histologyCorr %in% c("All",histo()),"")
      )
      if (input$dataset == "TCGA Lgg") {
        validate(need(input$colorBy != "Subtype" & input$separateBy != "Subtype", "Subtype available for GBM samples only")) 
      }
      myCorrTest(exprs(), input$gene, input$gene2, input$histologyCorr, input$subtype, separateBy = separateByInput())
    })
    
    #' Download the corrPlot
    output$downloadcorrPlot <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_corrPlot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file)
        myCorggPlot (exprs(), input$gene, input$gene2, input$histologyCorr, input$subtype, 
                     colorByInput(), separateByInput())
        dev.off()
      }
    )
    
    #' Generate an HTML table view of the data
    dataTable <- reactive({
      mRNA <- exprs()[ , input$gene, drop = FALSE]
      names(mRNA) <- paste0(input$gene,"_mRNA")
      data <- cbind(pDatas(),mRNA)
      if (input$dataset == "TCGA GBM" | input$dataset == "TCGA Lgg") {
        CN_status <- cnas()[,input$gene, drop = FALSE]
        CN_status[,1] <- factor(CN_status[,1], levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp")) 
        names(CN_status) <- paste0(input$gene,"_CN_status") 
        data <- cbind(data, CN_status)
      }
      data <- data[,colSums(is.na(data)) < nrow(data)] 
      data
    })
    
    #' Generate an HTML table view of the data
    output$table <- renderDataTable({
      if (input$gene == "")
        return(pDatas())
      dataTable()
    }, options = list(orderClasses = TRUE, lengthMenu = c(10, 30, 50), pageLength = 10))
    
    #' Download the data
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, ".csv", sep="")
      },
      content = function(file) {
        write.csv(dataTable(), file)
      }
    )
    
    #' Generate a graphic summary of the dataset, using rCharts
    output$piePlots <- renderUI({
      data <- pDatas()[,c("Histology","Grade","Recurrence","Subtype")]
      data <- data[,colSums(is.na(data)) < nrow(data)] # Removing unavailable (all NA) groups
      plot_output_list <- lapply(names(data), function(i) {
        plotname <- paste("plot", i, sep="")
        chartOutput(plotname, "nvd3") # "highcharts" for h plots
      })
      # Convert the list to a tagList - this is necessary for the list of items to display properly.
      do.call(tagList, plot_output_list)
    })
    
    observe ({                                                               
      data <- pDatas()[,c("Histology","Grade","Recurrence","Subtype")]
      data <- data[,colSums(is.na(data)) < nrow(data)]
      # Call renderChart for each one. 
      for (i in names(data)) {                                                    
        local({
          my_i <-i
          plotname <- paste("plot", my_i, sep="")
          output[[plotname]] <- renderChart2({
            plotData <- data.frame(table(data[, my_i]))
            n1 <- nPlot(x = "Var1", y = "Freq", data = plotData, type = "pieChart")
            n1$addParams(height = 400, width = 400)      
            #             n1$chart(showLegend = FALSE)
            n1
          })
        })
      }
    })
    
    #' Generate survival groups stratified by Histology, etc.
    output$survPlots <- renderUI({
      validate(
        need(input$dataset!= "Bao" & input$dataset!= "Reifenberger" & input$dataset!= "Gill", "Sorry, no survival data are available for this dataset")
      )
      df <- pDatas()[,1:7]
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      df <- droplevels.data.frame(subset(df, Histology!="Non-tumor")) # Exclude normal sample, not displaying properly
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      plot_output_list <- lapply(groups, function(i) {
        plot_surv_name <- paste("plotSurv", i, sep = "")
        plotOutput(plot_surv_name, height = 400, width = 400)
      })
      do.call(tagList, plot_output_list)
    })  
    
    observe({   
      df <- pDatas()[,1:7]
      df <- df[,colSums(is.na(df)) < nrow(df)] 
      df <- droplevels.data.frame(subset(df, Histology != "Non-tumor")) 
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      for (i in groups) {                                                    
        local({
          my_Survi <- i
          plot_surv_name <- paste("plotSurv", my_Survi, sep="")
          output[[plot_surv_name]] <- renderPlot({
            if (input$dataset== "Bao" || input$dataset == "Reifenberger" || input$dataset == "Gill") {
              return()
            }
            surv.status <- df[ ,"status"]
            surv.time <- df[ ,"survival"]
            my.Surv <- Surv(surv.time, surv.status == 1)
            expr.surv <- survfit(my.Surv ~ df[ ,my_Survi], data = df, conf.type = "none")
            plot(expr.surv, xlab = "Survival time (Months)", ylab = "% Surviving", 
                 yscale = 100,  col = 1:length((levels(df[ ,my_Survi]))), mark.time = FALSE, main = paste(my_Survi))
            legend("topright", legend = levels(df[ ,my_Survi]), col = 1:length((levels(df[ ,my_Survi]))), lty = 1)
          })
        })
      }
    })
    
    #' Reactive function to generate subtype call to pass to data table and download handler
    svm.call <- reactive ({
      # input$upFile will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      if (input$tumorType == "gbm") {
        train <- gbm.tcga[["expr"]]
      } else if (input$tumorType == "lgg") {
        train <- gorovets[["expr"]]
      }
      row.names(train) <- train[,"Sample"]
      train.exp <- data.frame(t(train[,-c(1:8)]))
      row.names(upData) <- upData[,"Sample"]
      upData <- upData [,-1]
      learn.exp <-t(upData)
      # Common genes of the two datasets
      genes <-intersect(row.names(train.exp),row.names(learn.exp))  
      # Subset the initial training matrix to only contain common genes
      train.1 <- train.exp[genes,]
      train.1 <- train.1[!is.na(train.1[,1]),]  # Check for NA's
      df.train <- train.1 - rowMeans(train.1)
      # Subset the learning matrix to only contain common genes
      learn.1 <- learn.exp[genes,]
      learn.1 <- learn.1[!is.na(learn.1[,1]),]
      df.learn <- learn.1- rowMeans(learn.1)
      Training <- as.factor(as.character(train$Subtype))
      require(kernlab)
      svm <- ksvm(t(df.train), #Training matrix
                  Training, # "Truth" factor
                  cross=10, # 10 fold cross validation to learn which SVM is best
                  kernel="vanilladot",  # vanilladot kernel keeps things simple
                  family="multinomial", # we are predicting 5 classes
                  prob.model=TRUE,  # We want a probability model included so we can give each predicted sample a score, not just a class
                  scale=FALSE)  # We already scaled before by mean-centering the data.
      subtype.call <- as.matrix(predict(svm, t(df.learn)))
      prob <- as.matrix(predict(svm, t(df.learn ), type="probabilities"))
      svm.call <- data.frame(subtype.call, prob)
      rownames(svm.call) <- rownames(t(df.learn))
      svm.call[,"Sample"] <- rownames(svm.call)
      svm.call
    })
    
    #' Rerndering the subtype call as a data table
    output$svm <- renderDataTable({ 
      if (is.null(input$upFile) | input$goSvm == 0)
        return(NULL)
      input$goSvm
      isolate({
        svm <- svm.call()
      })
    })
    
    #' Download the subtype call
    output$downloadSvm <- downloadHandler(
      filename = function() {
        paste("SVM_Subtype_Call.csv", sep="")
      },
      content = function(file) {
        write.csv(svm.call(), file)
      }
    )
    
    #' Reactive function to generate ssGSEA call to pass to data table and download handler
    gsva.call <- reactive ({
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      if (input$tumorType == "gbm") {
        gene_list <- subtype_list
      } else if (input$tumorType == "lgg") {
        return(NULL)
      }
      rownames(upData) <- upData$Sample
      exprs <- data.frame(t(upData[,-1]))
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq=FALSE,
                           min.sz=0, max.sz=10000, verbose=FALSE)
      subtype_scores <- t(gsva_results)
      subtype_final <- data.frame(subtype_scores, Subtype = names(subtype_list)[apply(subtype_scores,1,which.max)], 
                                  Sample = rownames(subtype_scores))
      subtype_final
    })
    
    #' Rerndering the subtype call as a data table
    output$gsva <- renderDataTable({ 
      if (is.null(input$upFile) | input$goGsva == 0)
        return(NULL)
      input$goGsva
      isolate({
      gsva <- gsva.call()
      })
    })
    
    #' Download the subtype call
    output$downloadGsva <- downloadHandler(
      filename = function() {
        paste("GSVA_Subtype_call.csv", sep="")
      },
      content = function(file) {
        write.csv(gvsa.call(), file)
      }
    )
    
    #' Generate the correlation table ##  corFast works locally but not on shinyapps.io
    corr <- reactive ({
      corr <- getCorr(datasetInputCor()[["expr"]], input$geneCor, input$histologyCorrTable)
      corr <- corr[order(-abs(corr$r)), ]
      corr
    })
    
    #' Generate a reactive element of the the correlation data 
    corrData <- reactive({  
          corr.table <- suppressWarnings(corr())  # suppressWarnings  is used to prevent the warning messages in the LGG dataset  
          if (input$sign == 0.01){
            corr.table <- subset(corr.table, p.value <= 0.01)
          } 
          if (input$cor == "Positive"){
            corr.table <- subset(corr.table, r > 0)
            corr.table <- corr.table[order(-corr.table$r), ]
          }
          if (input$cor == "Negative"){
            corr.table <- subset(corr.table, r < 0)
            corr.table <- corr.table[order(corr.table$r), ]
          } else {
            corr.table <- subset(corr.table, p.value <= 0.05)
          }
          corr.table
    })
    
    #' Generate an HTML table view of the correlation table 
    output$corrData <- renderDataTable({
      if (input$geneCor == "" | input$goCor == 0)
        return()
      validate(
        # Not all genes are available for all the dataset
        need(input$geneCor %in% names(exprs()),"Gene not available for this dataset")
      )  
      input$goCor
      isolate({  # https://groups.google.com/forum/#!searchin/shiny-discuss/submit$20button/shiny-discuss/3eXElZxZoaM/QtGCl-4qXzsJ      
      corr.table <- corrData()
      })
    }, options = list(orderClasses = TRUE))
    
    #' Download the correlation table 
    output$downloadCorrData <- downloadHandler(
      filename = function() {
        paste(input$geneCor, input$datasetCor, input$histologyCorrTable, "corrData.csv", sep="_")
      },
      content = function(file) {
        write.csv(corrData(),file)
      }
    )
    
  })
