library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)
library (shinysky)
library (shinyBS)
library (dplyr)

source("helpers.R")

#######################################
############## Datasets  ##############
#######################################
gbm.tcga <- readRDS("data/TCGA.GBM.Rds")
lgg.tcga <- readRDS("data/TCGA.LGG.Rds")
rembrandt <- readRDS("data/Rembrandt.Rds")
freije <- readRDS("data/Freije.Rds")
gravendeel <- readRDS("data/Gravendeel.Rds")
murat <- readRDS("data/Murat.Rds")
phillips <- readRDS("data/Phillips.Rds")
gene_names <- readRDS("data/gene_names.Rds")

#######################################
############## server.R  ##############
#######################################
shinyServer(
  function(input, output, session) {
    options(shiny.maxRequestSize=30*1024^2)    
    # Return the requested dataset
    datasetInput <- reactive({
      switch(input$dataset, 
             "TCGA GBM" = gbm.tcga,
             "TCGA Lgg" = lgg.tcga,
             "Rembrandt" = rembrandt,
             "Gravendeel" = gravendeel,
             "Phillips" = phillips,
             "Murat" = murat,
             "Freije" = freije)
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
             "Freije" = freije)
    })
    
    # Expression data
    exprs <- reactive({
      datasetInput()[["expr"]]
    })
    
    # Phenotype data
    pDatas <- reactive({
      datasetInput()[["pData"]]
    })
    
    # CNA data
    cnas <- reactive({
      datasetInput()[["cna"]]
    })
    
    # Text matching with the gene names list
    updateSelectizeInput(session, inputId = "gene", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "geneCor", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "gene2", choices = gene_names, server = TRUE) 

    # Return the requested plotType
    plotType <- reactive({
      switch(input$plotTypeSel, 
             Histology = "Histology", 
             "Copy number" = "Copy number", 
             Subtype = "Subtype", 
             Grade = "Grade", 
             Recurrence = "Recurrence")
    }) 
    
    # When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plotSelected <- reactive ({
      if (input$plotTypeSel %in% datasetInput()[["plotType"]]){ 
        plotType()
      } else {
        NULL
      }
    })
    
    # This will change the plot type available for a specific dataset
    observe({
    updateSelectInput(session, inputId = "plotTypeSel", choices = datasetInput()[["plotType"]], selected = plotSelected()) 
    }, priority = 10)
    
    # Caption with gene and dataset
    output$caption <- renderText({
      if (input$gene == "" )
        return()
      title <- paste(input$gene, "in", input$dataset, "dataset")
    })
    
#     # Help popup NOT WORKING YET
    output$help <- renderUI({ 
      helpPopup(title = "Help me pleaseeeeee", 
                content = includeMarkdown("tools/help.Rmd"), 
                placement = "right", 
                trigger = "click") 
    })
    
#     # Help popup alternative NOT WORKING YET
#         output$help <- renderUI({ 
#           helpModal(title = "Help me pleaseeeeee", link = "helpLink", content = includeMarkdown("tools/help.Rmd"))
#         })
    
    # Return the available histology, to be used in the updateSelectInput for correlation and survival
    histo <- reactive({
      levels(datasetInput()[["pData"]][,"Histology"])
    })

    histoCor <- reactive({
      levels(datasetInputCor()[["pData"]][,"Histology"])
    })

    observe({
      # This will change the value of input$histologyCorrTable, based on histological group available for that dataset
      updateSelectInput(session, inputId = "histologyCorrTable", choices = c("All", histoCor()), selected = "All")
    })
    
    # When switching datasets, if the selected histo is not available it will choose GBM (the last histo of the list)
    histoSurvSelected <- reactive ({
      if (input$histologySurv %in% histo()){
        input$histologySurv
      } else {
        tail(histo(), n=1)
      }
    })

    # When switching datasets, if the selected histo is not available it will choose "All"
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
    
    
    # Need a wrapper around the densityClick input so we can manage whether or 
    # not the click occured on the current Gene. If it occured on a previous
    # gene, we'll want to mark that click as 'stale' so we don't try to use it
    # later. https://gist.github.com/trestletech/5929598
    currentClick <- list(click=NULL, stale=FALSE)

    handleClick <- observe({
      if (!is.null(input$densityClick) && !is.null(input$densityClick$x)){
        currentClick$click <<- input$densityClick
        currentClick$stale <<- FALSE
      }
    }, priority=100)
    
    getCutoff <- reactive({
      # We need this function to subscribe to a couple of dependencies. Without
      # explicitly providing these two at the top, this function may return a 
      # cached value without realizing that a new click has occured, or a new 
      # gene has been loaded.
      input$densityClick
      geneExp()        
      # See if there's been a click since the last gene change.
      if (!is.null(currentClick$click) && !currentClick$stale){
        return(currentClick$click$x)
      }       
      return (median(geneExp()))
    })
    
    # Extract the relevant GBM expression values.
    geneExp <- reactive({
      df <- subset(exprs(), Histology == "GBM")
      if (input$gcimp){
        df <- subset(exprs(), Subtype != "G-CIMP")
      }
      geneExp <- df[ ,input$gene]
      currentClick$stale <<- TRUE
      geneExp
    })
    
    # Extract the Hazard ratios for the input gene.
    HR <- reactive ({
      HR <- getHR(exprs(), input$gene, input$gcimp)
    })
    
    #' Render a plot to show the the Hazard ratio for the gene's expression values
    output$hazardPlot <- renderPlot({        
      if (input$gene == "")
        return()
      validate(
        need(histoSurvSelected()  == "GBM", "Interactive HR plot currently available only for GBM samples")
      )
      # Wrap the entire expensive operation with withProgress 
      withProgress(session, min = 1, max = 5, {
        setProgress(message = "Calculating, please wait",
                    detail = "This may take a few moments...")
        for (i in 1:5) {
          setProgress(value = i)
          Sys.sleep(0.5)
        }          
        # Plot the hazardplot 
        hazardPlot(HR(), input$quantile)
        # Add a vertical line to show where the current cutoff is.
        abline(v = getCutoff(), col = 4)
      })
    }, bg = "transparent")
    
    # A reactive survival formula
    survivalFml <- reactive({
      df <- subset (exprs(), Histology == "GBM")
      if (input$gcimp){
        df <- subset (exprs(), Subtype != "G-CIMP")
      }
      # Create the groups based on which samples are above/below the cutoff
      expressionGrp <- as.integer(geneExp() < getCutoff())
      # Create the survival object
      surv <- with(df, Surv(survival, status == 1))
      return(surv ~ expressionGrp)
    })
    
    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      if (input$gene == "")
        return()
      validate(
        need(histoSurvSelected()  == "GBM","Interactive HR plot currently available only for GBM samples")
      )
      cutoff <- getCutoff()
      surv <- survivalFml()
      kmPlot(cutoff, surv)
    })

    #' Create a Kaplan Meier plot with cutoff based on quantiles
    output$survPlot <- renderPlot({
      if (input$gene == "" | input$histologySurv == "")
        return()
      validate(
        need(input$histologySurv %in% histo(),"")
      )# Trying to avoid an error when switching datasets in case the choosen histology is not available.
      validate(
        need(input$histologySurv != "Non-tumor","Sorry, no survival data are available for this group")
      )
      survivalPlot (exprs(), input$gene, group = input$histologySurv, cutoff = input$cutoff, subtype = input$subtypeSurv, 
                    gcimp = input$gcimpSurv)
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({     
      if (input$gene == "" )
        return()
      validate(
        need(plotType() %in% datasetInput()[["plotType"]],"")
      ) # Trying to avoid an error when switching datasets in case the plotType is not available.
      validate(
        need(input$gene %in% names(exprs()),"Gene not available for this dataset")
      ) # Not all genes are available for all the dataset
      ggboxPlot(exprs = exprs(), cna = cnas(), gene = input$gene, plotType = plotType(), scale = input$scale, 
                stat = input$stat, colBox = input$colBox, colStrip = input$colStrip, bw = input$bw) 
      # I needed to create the ggboxPlot function (see helper file) to use it in the output$downloadPlot .... OTHER WAY TO DO IT??
    })
    
    #' Tukey post-hoc test
    output$tukeyTest <- renderPrint(width = 800, {    
      if (input$gene == "" )
        return()
      validate(
        need(plotType() %in% datasetInput()[["plotType"]],"")
      ) # Trying to avoid an error when switching datasets in case the plotType is not available.
      mRNA <- exprs()[ ,input$gene]
      if (plotType() == "Copy number") {
        group <- cnas()[ ,input$gene]
        group <- factor(group, levels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        group <- droplevels(group)
      } else {
        group <- exprs()[ ,plotType()]
      }
      if (any(!is.na(group))) {
        data <- data.frame(mRNA, group)
        data <- na.omit(data)
      }
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff), ]
      tukey
     })
    
    
    #' Pairwise t test
    output$pairwiseTtest <- renderPrint({     
      if (input$gene == "" )
        return()
      validate(
        need(plotType() %in% datasetInput()[["plotType"]],"")
      ) # Trying to avoid an error when switching datasets in case the plotType is not available.
      mRNA <- exprs()[ ,input$gene]
      if (plotType() == "Copy number") {
        group <- cnas()[ ,input$gene]
        group <- factor(group, levels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        group <- droplevels(group)
      } else {
        group <- exprs()[ ,plotType()]
      }
      if (any(!is.na(group))) {
        data <- data.frame(mRNA, group)
        data <- na.omit(data)
      }
      pttest <- pairwise.t.test(mRNA, group, na.rm= TRUE, p.adj = "bonferroni", paired = FALSE)[[3]]
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
    
    # Get the download dimensions.
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
      # The argument content below takes filename as a function
      # and returns what's printed to it.
      content = function(con) {
        # Gets the name of the function to use from the 
        # downloadFileType reactive element. Example:
        # returns function pdf() if downloadFileType == "pdf".
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(con, width = downloadPlotWidth(), height = downloadPlotHeight())
        ggboxPlot(exprs = exprs(), cna = cnas(), gene = input$gene, plotType = plotType(), scale = input$scale, 
                  stat = input$stat, colBox = input$colBox, colStrip = input$colStrip, bw = input$bw)
        dev.off(which=dev.cur())
      }
    )
    
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
    
    #' Generate the correlation table ##  CURRENTLY TOO SLOW
    corr <- reactive ({
      corr <- getCorr(datasetInputCor()[["expr"]], input$geneCor, input$histologyCorrTable)
      corr <- corr[order(-abs(corr$r)), ]
      corr
    })

    #' Generate an HTML table view of the the correlation table 
    output$corrData <- renderDataTable({
      if (input$geneCor == "" | input$goCor == 0)
        return()
      isolate({  # https://groups.google.com/forum/#!searchin/shiny-discuss/submit$20button/shiny-discuss/3eXElZxZoaM/QtGCl-4qXzsJ
        withProgress(session, min=1, max=5, {
          setProgress(message = "Calculating, please wait",
                      detail = "Be patient, this takes forever...")
          for (i in 1:5) {
            setProgress(value = i)
            Sys.sleep(0.5)
          }          
          corr.table <- suppressWarnings(corr())  # suppressWarnings  is used to prevent the warning messages in the LGG dataset  
          if (input$sign == 0.01){
            corr.table <- subset(corr.table, p <= 0.01)
          } 
          if (input$cor == "Positive"){
            corr.table <- subset(corr.table, r > 0)
            corr.table <- corr.table[order(-corr.table$r), ]
          }
          if (input$cor == "Negative"){
            corr.table <- subset(corr.table, r < 0)
            corr.table <- corr.table[order(corr.table$r), ]
          } else {
            corr.table <- subset(corr.table, p <= 0.05)
          }
          corr.table
        })
      })
    })

    #' Download the correlation table 
    output$downloadCorrData <- downloadHandler(
      filename = function() {
        paste(input$geneCor, input$datasetCor, input$histologyCorrTable, "corrData.csv", sep="_")
      },
      content = function(file) {
        write.csv(corr(),file)
      }
    )
    
    ## Generate reactive Inputs for the corrPlot to be used also to download the complete plot
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
      if (input$gene == "")
        return()
      validate(
        need(input$gene2 != "", "Please enter Gene 2")
      )
      validate(
        need(input$histologyCorr %in% c("All",histo()),"")
      ) # Trying to avoid an error when switching datasets in case the choosen histology is not available.
      myCorggPlot(exprs(), input$gene, input$gene2, input$histologyCorr, input$subtype, 
                   colorBy = colorByInput(), separateBy = separateByInput())
    })
    
    #' Generate a summary of the correlation test
    output$corrTest <-renderPrint({     
      if (input$gene == "" | input$gene2 == "")
        return()
      validate(
        need(input$histologyCorr %in% c("All",histo()),"")
      ) # Trying to avoid an error when switching datasets in case the choosen histology is not available.
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

    #' Generate a graphic summary of the dataset, using rCharts
    output$piePlots <- renderUI({
      data <- pDatas()[,c("Histology","Grade","Recurrence","Subtype")]
      data <- data[,colSums(is.na(data)) < nrow(data)]
      plot_output_list <- lapply(names(data), function(i) {
        plotname <- paste("plot", i, sep="")
        chartOutput(plotname, "highcharts")
      })
      # Convert the list to a tagList - this is necessary for the list of items to display properly.
      do.call(tagList, plot_output_list)
    })
    
    observe ({                                                               
      data <- pDatas()[,c("Histology","Grade","Recurrence","Subtype")]
      data <- data[,colSums(is.na(data)) < nrow(data)]
      # Call renderChart for each one. 
      for (i in names(data)) {                                                    
        # Need local so that each item gets its own name. Without it, the value
        # of i in the renderChart2() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <-i
          plotname <- paste("plot", my_i, sep="")
          output[[plotname]] <- renderChart2({
            plotData <- data.frame(table(data[, my_i]))
            hPlot(x = "Var1", y = "Freq", data = plotData, type = "pie", title = my_i)
          })
        })
      }
    })

    #' Generate an HTML table view of the data
    output$table <- renderDataTable({
      if (input$gene == "")
        return(pDatas())
      data <- getData(exprs(), input$gene)
      data.frame(data)
    })
    
    #' Download the data
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, ".csv", sep="")
      },
      content = function(file) {
        write.csv(getData(exprs(), input$gene), file)
      }
    )

    #' Generate survival groups stratified by Histology, etc.
    output$survPlots <- renderUI({
      df <- pDatas()
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      df <- droplevels.data.frame(subset(df, Histology!="Non-tumor")) # Exclude normal sample, not displaying properly
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      plot_output_list <- lapply(groups, function(i) {
        plot_surv_name <- paste("plotSurv", i, sep = "")
        plotOutput(plot_surv_name, height = 400, width = 400)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
    })  
    
    observe({   
      df <- pDatas()
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      df <- droplevels.data.frame(subset(df, Histology!="Non-tumor")) # Exclude normal sample, not displaying properly
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      for (i in groups) {                                                    
        # Need local so that each item gets its own name. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_Survi <- i
          plot_surv_name <- paste("plotSurv", my_Survi, sep="")
          output[[plot_surv_name]] <- renderPlot({
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
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                         quote=input$quote)
      if (input$svm == "gbm") {
      tcga <- gbm.tcga[["expr"]]
      row.names(tcga) <- tcga[,"Sample"]
      TCGA.exp <- data.frame(t(tcga[,-c(1:8)]))
      row.names(upData) <- upData[,"Sample"]
      upData <- upData [,-1]
      df.exp <-t(upData)
      genes <-intersect(row.names(TCGA.exp),row.names(df.exp))  # common genes of the two datasets
      ## Subset the initial TCGA matrix to only contain common genes
      TCGA.1 <- TCGA.exp[genes,]
      TCGA.1 <- TCGA.1[!is.na(TCGA.1[,1]),]  #Check for NA's
      TCGA.train <- TCGA.1 - rowMeans(TCGA.1)
      ## Subset the inputDf matrix to only contain common genes
      df.1 <- df.exp[genes,]
      df.1 <- df.1[!is.na(df.1[,1]),]
      df.learn <- df.1- rowMeans(df.1)
      Training <- as.factor(as.character(tcga$Subtype))
      require(kernlab)
      svm.TCGA <- ksvm(t(TCGA.train), #Training matrix
                       Training, # "Truth" factor
                       cross=10, # 10 fold cross validation to learn which SVM is best
                       kernel="vanilladot",  # vanilladot kernel keeps things simple
                       family="multinomial", # we are predicting 5 classes
                       prob.model=TRUE,  # We want a probability model included so we can give each predicted sample a score, not just a class
                       scale=FALSE)  # We already scaled before by mean-centering the data.
      subtype.call      <- as.matrix(predict(svm.TCGA, t(df.learn)))
      prob      <- as.matrix(predict(svm.TCGA, t(df.learn ), type="probabilities"))
      svm.call       <- data.frame(subtype.call ,prob)
      rownames(svm.call) <- rownames(t(df.learn))
      svm.call[,"Sample"] <- rownames(svm.call)
      svm.call}
      else if (input$svm == "lgg") {
        return(NULL) # LGG subtype not yet active
      }
      svm.call
    })

    #' Rerndering the subtype call as a data table
    output$svm <- renderDataTable({ 
      inFile <- input$upFile
      if (is.null(inFile))
        return(NULL)
      # Wrap the entire expensive operation with withProgress 
      withProgress(session, min = 1, max = 5, {
        setProgress(message = "Calculating, please wait",
                    detail = "Be patients, switching to another tab will crash GlioVis ...")
        for (i in 1:5) {
          setProgress(value = i)
          Sys.sleep(0.5)
        }
        if (input$svm == "gbm") {
          svm <- svm.call()
        } else if (input$svm == "lgg") {
          return(NULL) # LGG subtype not yet active
        }
        svm
      })
    })

    #' Download the subtype call
    output$downloadSvm <- downloadHandler(
      filename = function() {
        paste("My_SubtypeCall.csv", sep="")
      },
      content = function(file) {
        write.csv(svm.call(), file)
      }
    )

})
