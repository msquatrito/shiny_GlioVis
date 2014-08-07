library (shiny)
library (survival)
library (weights)
library (shinyIncubator)
library (ggplot2)
library (gridExtra)
library (rCharts)

source("helpers.R")
gbm.tcga <- readRDS("data/TCGA.GBM.Rds")
lgg.tcga <- readRDS("data/TCGA.LGG.Rds")
rembrandt <- readRDS("data/Rembrandt.Rds")
freije <- readRDS("data/Freije.Rds")
gravendeel <- readRDS("data/Gravendeel.Rds")
murat <- readRDS("data/Murat.Rds")
phillips <- readRDS("data/Phillips.Rds")
gene_names <- readRDS("data/gene_names.Rds")


shinyServer(
  function(input, output, session) {
    
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
    
    # Return the requested plotType
    plotType <- reactive({
      switch(input$plotTypeSel, 
             Histology = "Histology", 
             "Copy number" = "Copy number", 
             Subtype = "Subtype", 
             Grade = "Grade", 
             Recurrence = "Type")
    })
    
    # Text matching with the gene names list
    updateSelectizeInput(session, "gene",  choices = gene_names, server = TRUE)
    
    # Caption with gene and dataset
    output$caption <- renderText({
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      title <- paste(input$gene, "in", input$dataset, "dataset")
    })
    
#     addTooltip(session, "dataset", "Select a dataset", "right", trigger="hover")
#     addTooltip(session, "gene", "Enter gene name", "right", trigger="hover")
#     addTooltip(session, "plotTypeSel", "Select one of the available plot", "right", trigger="hover")
    
    # Help popup NOT WORKING YET
    output$help <- renderUI({ 
      helpPopup(title = "Help me pleaseeeeee", 
                #                 content = "", 
                content = includeMarkdown("tools/help.Rmd"), # I've still to create the Rmd file
                placement = "right", trigger = "click") 
    })
    
    # Help popup alternative NOT WORKING YET
    #     output$help <- renderUI({ 
    #       helpModal("Help me pleaseeeeee","help",includeMarkdown("tools/help.Rmd"))
    #     })
    #     
    
    # Return the available histology, to be used in the updateSelectInput for correlation and survival
    histo <- reactive({
      levels(datasetInput()[["pData"]][,"Histology"])
    })
    
    observe({
      # This will change the value of input$gene1, based on input$gene
      updateSelectizeInput(session, "gene1", choices = gene_names, server = TRUE, selected = input$gene)
      # This will change the value of input$gene2, based on input$gene
      updateSelectizeInput(session, "gene2", choices = gene_names, server = TRUE, selected = input$gene)
      # This will change the value of input$histologySurv, based on histological group available for that dataset
      updateSelectInput(session, "histologySurv", choices = histo(), selected = tail(histo(), n=1))
      # This will change the value of input$histologyCorr, based on histological group available for that dataset
      updateSelectInput(session, "histologyCorr", choices = c("All", histo()), selected = "All")
      # This will change the value of input$histologyCorr, based on histological group available for that dataset
      updateSelectInput(session, "histologyCorrTable", choices = c("All", histo()), selected = "All")
      # This will change the plot type available for a specific dataset
      updateSelectInput(session, "plotTypeSel", choices = datasetInput()[["plotType"]], selected = NULL) # selected = datasetInput()[["plotType"]][1]
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
      df <- datasetInput()[["expr"]]
      df <- subset (df, Histology == "GBM")
      if (input$gcimp){
        df <- subset (df, Subtype != "G-CIMP")
      }
      geneExp <- df[ ,input$gene]
      currentClick$stale <<- TRUE
      geneExp
    })
    
    # Extract the Hazard ratios for the input gene.
    HR <- reactive ({
      df <- datasetInput()[["expr"]]
      HR <- getHR(df, input$gene, input$gcimp)
    })
    
    #' Render a plot to show the distribution of the gene's expression
    output$hazardPlot <- renderPlot({        
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      # Wrap the entire expensive operation with withProgress 
      withProgress(session, min=1, max=5, {
        setProgress(message = "Calculating, please wait",
                    detail = "This may take a few moments...")
        for (i in 1:5) {
          setProgress(value = i)
          Sys.sleep(0.5)
        }          
        # Plot the hazardplot
        HR <- HR()  
        hazardPlot(HR, input$quantile)
        # Add a vertical line to show where the current cutoff is.
        abline(v=getCutoff(), col=4)
        }
      )
    }, bg="transparent")
    
    # A reactive survival formula
    survivalFml <- reactive({
      df <- datasetInput()[["expr"]]
      df <- subset (df, Histology == "GBM")
      if (input$gcimp){
        df <- subset (df, Subtype != "G-CIMP")
      }
      # Create the groups based on which samples are above/below the cutoff
      expressionGrp <- as.integer(geneExp() < getCutoff())
      # Create the survival object
      surv <-with(df, Surv(survival,status== 1))
      return(surv ~ expressionGrp)
    })
    
    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      cutoff <- getCutoff()
      surv <- survivalFml()
      kmPlot(cutoff,surv)
    })

    #' Create a Kaplan Meier plot with cutoff based on quantiles
    output$survPlot <- renderPlot({
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" | input$histologySurv == "" )
        return()
      df <- datasetInput()[["expr"]]
      survivalPlot (df, input$gene, group = input$histologySurv, cutoff = input$cutoff, subtype = input$subtypeSurv, gcimp = input$gcimpSurv)
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({     
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      exprs <- datasetInput()[["expr"]]
      cna <- datasetInput()[["cna"]]
      
      ggboxPlot(exprs = exprs, cna = cna, gene = input$gene, plotType = plotType(), scale = input$scale, 
                stat = input$stat, colBox = input$colBox, colStrip = input$colStrip) 
      # I needed to create the ggboxPlot function (see helper file) to use it in the output$downloadPlot .... OTHER WAY TO DO IT??
    })
    
    #' Tukey post-hoc test
    output$tukeyTest <- renderPrint({     
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      exprs <- datasetInput()[["expr"]]
      cna <- datasetInput()[["cna"]]
      mRNA <- exprs[ ,input$gene]
      if (plotType() == "Copy number") {
        group <- cna[, input$gene]
        group <- factor(group, levels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        group <- droplevels(group)
      } else {
        group <- exprs[ ,plotType()]
      }
      if (any(!is.na(group))) {
        data <- data.frame(mRNA,group)
        data <- na.omit(data)
      }
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj,p.levels=c(.001, .01, .05, 1), 
                                                symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff),]
      tukey
     })
    
    
    #' pairwise t test
    output$pairwiseTtest <- renderPrint({     
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      exprs <- datasetInput()[["expr"]]
      cna <- datasetInput()[["cna"]]
      mRNA <- exprs[ ,input$gene]
      if (plotType() == "Copy number") {
        group <- cna[, input$gene]
        group <- factor(group, levels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        group <- droplevels(group)
      } else {
        group <- exprs[ ,plotType()]
      }
      if (any(!is.na(group))) {
        data <- data.frame(mRNA,group)
        data <- na.omit(data)
      }
      pttest <- pairwise.t.test(mRNA, group, p.adj="bonferroni", paired=F)[[3]]
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
        ggboxPlot(exprs = datasetInput()[["expr"]], cna = datasetInput()[["cna"]], gene = input$gene, 
                  plotType = plotType(), scale = input$scale, stat = input$stat, colBox = input$colBox, colStrip = input$colStrip)
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
        survivalPlot (datasetInput()[["expr"]], input$gene, group = input$histologySurv, 
                      cutoff = input$cutoff, subtype = input$subtypeSurv, gcimp = input$gcimpSurv)
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
      corr <- getCorr(datasetInput()[["expr"]], input$gene, input$histologyCorrTable)
      corr <- corr[order(-abs(corr$r)), ]
      corr
    })

    #' Generate an HTML table view of the the correlation table 
    output$corrData <- renderDataTable({
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return()
      # Wrap the entire expensive operation with withProgress 
      withProgress(session, min=1, max=5, {
        setProgress(message = "Calculating, please wait",
                    detail = "This takes forever...")
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

    #' Download the correlation table 
    output$downloadCorrData <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_corrData.csv", sep="")
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
      if (input$gene1 == "Enter gene, eg: EGFR" | input$gene1 == "" |
            input$gene2 == "Enter gene, eg: EGFR" | input$gene2 == "")
        return()
      df <- datasetInput()[["expr"]]
      myCorggPlot (df, input$gene1, input$gene2, 
                   input$histologyCorr, input$subtype , 
                   colorBy = colorByInput (), separateBy = separateByInput())
    })
    
    #' Generate a summary of the correlation test
    output$corrTest <-renderPrint({     
      if (input$gene1 == "Enter gene, eg: EGFR" | input$gene1 == "" |
            input$gene2 == "Enter gene, eg: EGFR" | input$gene2 == "")
        return()
      df <- datasetInput()[["expr"]]  
      myCorrTest (df, input$gene1, input$gene2, input$histologyCorr, input$subtype)
    })
    
    #' Download the corrPlot
    output$downloadcorrPlot <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, "_corrPlot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file)
        myCorggPlot (datasetInput()[["expr"]], input$gene1, input$gene2, 
                     input$histologyCorr, input$subtype, 
                     colorByInput (), separateByInput())
        dev.off()
      }
    )
    
#     #' Generate a summary of the dataset NOT SURE IS USEFUL
#     output$summary <- renderPrint({
#       if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
#         return()
#       data <- getData (datasetInput()[["expr"]], input$gene)
#       summary(data[,-c(1,7)])
#     })

    #' Generate a graphic summary of the dataset, using rCharts
    output$piePlots <- renderUI({
      data <- datasetInput()[["expr"]][,c("Histology","Grade","Type","Subtype")]
      data <- data[,colSums(is.na(data)) < nrow(data)]
      plot_output_list <- lapply(names(data), function(i) {
        plotname <- paste("plot", i, sep="")
        chartOutput(plotname, "highcharts")
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
    })
    
    observe ({                                                               
      data <- datasetInput()[["expr"]][,c("Histology","Grade","Type","Subtype")]
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
      if (input$gene == "Enter gene, eg: EGFR" | input$gene == "" )
        return(datasetInput()[["pData"]])
      data <- getData (datasetInput()[["expr"]], input$gene)
      data.frame(data)
    })
    
    #' Download the data
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$gene, "_", input$dataset, ".csv", sep="")
      },
      content = function(file) {
        write.csv(getData (datasetInput()[["expr"]], input$gene),file)
      }
    )

    #' Generate survival groups stratified by Histology, etc.
    output$survPlots <- renderUI({
      df <- datasetInput()[["pData"]]
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      df <- droplevels.data.frame(subset(df, Histology!="Non-tumor")) # Exclude normal sample, not displayng properly
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      plot_output_list <- lapply(groups, function(i) {
        survPlotname <- paste("plotSurv", i, sep="")
        plotOutput(survPlotname, height = 400, width = 400)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
    })  
    
    observe ({   
      df <- datasetInput()[["pData"]]
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      df <- droplevels.data.frame(subset(df, Histology!="Non-tumor")) # Exclude normal sample, not displayng properly
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      for (i in groups) {                                                    
        # Need local so that each item gets its own name. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_Survi <-i
          survPlotname <- paste("plotSurv", my_Survi, sep="")
          output[[survPlotname]] <- renderPlot({
            surv.status <- df[ ,"status"]
            surv.time <- df[ ,"survival"]
            my.Surv <- Surv(surv.time, surv.status == 1)
            expr.surv <- survfit(my.Surv ~ df[,my_Survi], data = df, conf.type = "none")
            plot(expr.surv, xlab = "Survival time (Months)", ylab = "% Surviving", 
                 yscale = 100,  col = 1:length((levels(df[,my_Survi]))), mark.time = FALSE,
                 main = paste(my_Survi))
            legend("topright", legend = levels(df[,my_Survi]), col = 1:length((levels(df[,my_Survi]))), lty = 1)
          })
        })
      }
    })


})
