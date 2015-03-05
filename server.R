# datasets
# gbm.tcga <- readRDS("data/TCGA.GBM.Rds")
# lgg.tcga <- readRDS("data/TCGA.LGG.Rds")
# rembrandt <- readRDS("data/Rembrandt.Rds")
# freije <- readRDS("data/Freije.Rds")
# gravendeel <- readRDS("data/Gravendeel.Rds")
# murat <- readRDS("data/Murat.Rds")
# phillips <- readRDS("data/Phillips.Rds")
# reifenberger <- readRDS("data/Reifenberger.Rds")
# bao <- readRDS("data/Bao.Rds")
# gill <- readRDS("data/Gill.Rds")
# gorovets <- readRDS("data/Gorovets.Rds")
# nutt <- readRDS("data/Nutt.Rds")
# ducray <- readRDS("data/Ducray.Rds")
# grzmil <- readRDS("data/Grzmil.Rds")
# donson <- readRDS("data/Donson.Rds")
# li <- readRDS("data/Li.Rds")
# vital <- readRDS("data/Vital.Rds")
# joo <- readRDS("data/Joo.Rds")
# oh <- readRDS("data/Oh.Rds")

# server.R for Gliovis
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
             "Nutt" = nutt,
             "Ducray" = ducray,
             "Grzmil" = grzmil,
             "Donson" = donson,
             "Li" = li,
             "Vital" = vital,
             "Joo" = joo,
             "Oh" = oh)
    })
    
    #' Switch the dataset for the correlation
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
             "Nutt" = nutt,
             "Ducray" = ducray,
             "Grzmil" = grzmil,
             "Donson" = donson,
             "Li" = li,
             "Vital" = vital,
             "Joo" = joo,
             "Oh" = oh)
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
    updateSelectizeInput(session, inputId = "genelist", choices = gene_names, server = TRUE) 
    
    #' Required for the conditional panel 'geneslist' to work correctly
    observe({
      if(input$tab1 != 3)
        updateTabsetPanel(session, inputId = "tabCorr", selected = "2genes")
    })
    
    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plotSelected <- reactive ({
      if (input$plotTypeSel %in% plotList[[input$dataset]]){ 
        input$plotTypeSel
      } else {
        NULL
      }
    })
    
    #' Return the names of the available user-defined plots
    plotUserSelection <- reactive ({
      data <- pDatas()[,!names(pDatas())%in%c("Sample","Histology","Grade","Recurrence","Subtype", "CIMP_status", "survival",
                                              "status", "Age", "ID","Patient_ID","Sample_ID", "Matching.sample", "Therapy_Class","title")] # Exlude pre-defined plots and numeric variables
      n <- colnames(data)
      n
    })
    
    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plotUserSelected <- reactive ({
      if (input$plotTypeUserSel %in% plotUserSelection()){ 
        input$plotTypeUserSel
      } else {
        NULL
      }
    })
    
    
    #' Change the plot type available for a specific dataset
    observe({
      updateSelectInput(session, inputId = "plotTypeSel", choices = plotList[[input$dataset]], selected = plotSelected()) 
    }, priority = 10)
    
    observe({
      updateSelectInput(session, inputId = "plotTypeUserSel", choices = plotUserSelection(), selected = plotUserSelected())
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
      if (input$histologySurv %in% c("All", histo())){
        input$histologySurv
      } else {
        "All"
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
      updateSelectInput(session, inputId = "histologySurv", choices = c("All", histo()), selected = histoSurvSelected())
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
    
    #' Generate a dataframe with the data to plot 
    data <- reactive({     
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the plotType is not available.
        need(plotSelected() %in% plotList[[input$dataset]],""),
        need(plotUserSelected() %in% plotUserSelection(),"")
      ) 
      mRNA <- exprs()[ ,input$gene]
      if (input$scale) {
        mRNA <- scale(mRNA)
      }
      data <- cbind(mRNA, exprs()[,2:6])
      data <- cbind(data, pDatas()[,!names(pDatas()) %in% names(data)])
      if (input$dataset == "TCGA GBM" | input$dataset == "TCGA Lgg") {
        Copy_number <- cnas()[ ,input$gene]
        Copy_number <- factor(Copy_number, levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp")) 
        data <- cbind(Copy_number,data)
      }
      if (input$primary & any(!is.na(data$Recurrence))) {
        data <- subset (data, Recurrence == "Primary")
      }
      data <- data[,colSums(is.na(data)) < nrow(data)] 
      data$group <- data[,plotType()]
      data <- subset(data,!is.na(group))
    })

    #' Generate radiobuttons for the various categories in the pData
    output$colorPoints <- renderUI({
      validate(need(input$gene != "",""))
      colnames <- names(data())[!names(data()) %in% c("group","mRNA","Sample","status","survival")]
      # Create the radiobuttons for the different pData categories
      radioButtons("colorP", "Color by:", choices  = colnames, selected = plotType())
    })
    
    plotType <- reactive({
      if (input$plotType == "Pre-defined"){
        plotSelected()
      } else if (input$plotType == "User-defined") {
        plotUserSelected()
      }
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({
    # To avoid an error when switching datasets in case the colStrip is not available.
      if(input$colStrip){
        colnames <- names(data())[!names(data()) %in% c("group","mRNA","Sample","status","survival")]
        validate(need(input$colorP %in% colnames,""))
      }
      ggboxPlot(data = data (), xlabel = plotType(), scale = input$scale, stat = input$tukeyPlot, colBox = input$colBox, 
                colStrip = input$colStrip, colorPoints = input$colorP, bw = input$bw) 
    })
       
    #' Summary statistic
    output$summary <- renderTable({    
      stat <- data.frame(data() %>%
                           group_by(group) %>%
                           summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data())[1],2), "%)" ), # prop.table
                                     median = median (mRNA, na.rm=T), mad = mad(mRNA, na.rm=T),mean = mean(mRNA, na.rm=T), sd = sd(mRNA, na.rm=T)))
      row.names(stat) <- stat$group
      stat <- stat[,-1]
      stat      
    }, align='rrrrrr')
    
    #' Tukey post-hoc test
    output$tukeyTest <- renderTable({    
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data()))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff), ]
      tukey
    })
    
    
    #' Pairwise t test
    output$pairwiseTtest <- renderTable({     
      data <- data()
      x <- data$mRNA # mRNA
      y <- data$group # group
      pttest <- pairwise.t.test(x, y, na.rm= TRUE, p.adj = "bonferroni", paired = FALSE)[[3]]
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
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_", input$plotTypeSel,  
               ".", downloadPlotFileType() )
      },     
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(file, width = downloadPlotWidth(), height = downloadPlotHeight())
        ggboxPlot(data = data (), xlabel = plotType(), scale = input$scale, stat = input$tukeyPlot, colBox = input$colBox, 
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
    
    #' Reactive expressions for the conditonal panels to work in the right way
    gcimpSurv <- reactive ({
      if(input$histologySurv == "GBM") {
        gcimpSurv <- input$gcimpSurv
      } else {
        gcimpSurv <- FALSE
      }
      gcimpSurv
    })
    
    primarySurv <- reactive ({
      if(input$histologySurv == "GBM") {
        primarySurv <- input$primarySurv
      } else {
        primarySurv <- FALSE
      }
      primarySurv
    })
    
    allSubSurv <- reactive ({
      if(input$histologySurv == "GBM" & input$subtypeSurv == "All") {
        allSubSurv <- input$allSubSurv
      } else {
        allSubSurv <- FALSE
      }
      allSubSurv
    })
    
    #' Extract the GBM samples for the interactive HR plot.
    survData <- reactive({    
      df <- exprs()
      df <- subset(df, !is.na(df$status))
      # subset to primary GBM, in case both primary and recurrent samples are available
      if (primarySurv() & any(!is.na(df$Recurrence))) {
        df <- subset (df, Histology == "GBM" & Recurrence == "Primary")
      }
      df <- subset (df, Histology == "GBM")
      if (input$subtypeSurv != "All") {
        df <- subset (df, Subtype == input$subtypeSurv)
      }
      # exclude G-CIMP is selected
      if (gcimpSurv()){
        df <- subset (df, CIMP_status != "G-CIMP")
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
        need(!input$dataset %in% noSurvDataset, "Sorry, no survival data are available for this dataset")%then%
          need(input$histologySurv != "Non-tumor","Sorry, no survival data are available for this group")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")
      )
    })
    
    #' busy indicator when switching surv tab
    #' http://stackoverflow.com/questions/18237987/show-that-shiny-is-busy-or-loading-when-changing-tab-panels
    output$activeTabSurv <- reactive({
      return(input$tabSurv)
    })
    outputOptions(output, 'activeTabSurv', suspendWhenHidden=FALSE)
    
    #' Render a plot to show the the Hazard ratio for the gene's expression values
    output$hazardPlot <- renderPlot({        
      validate(need(histoSurvSelected() == "GBM", "Interactive HR plot currently available only for GBM samples"))
      survNeed()
      validate(need(!input$dataset %in% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      input$tabSurv
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
      survNeed()
      validate(need(input$histologySurv %in% c("All", histo()),""))      
      df <- exprs()
      df <- subset(df, !is.na(df$status))
      if (input$histologySurv != "All"){
        df <- subset(exprs(), Histology == input$histologySurv)
      }
      if (input$histologySurv == "GBM" & input$subtypeSurv != "All") {
        df <- subset (df, Subtype == input$subtypeSurv)
      }
      if (gcimpSurv()){
        df <- subset (df, CIMP_status != "G-CIMP")
      }
      if (primarySurv() & any(!is.na(df$Recurrence))) {
        df <- subset (df, Recurrence == "Primary")
      }
      mRNA <- df[ ,input$gene]
      mRNA.values <- round(mRNA[!is.na(mRNA)],2)
      # Generate a vector of continuos values, excluding the first an last value
      mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) 
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxmRNA <- renderPlot({    
      validate(need(input$mInput, ""))      
      mRNA <- round(mRNAsurv(),2)
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|")
      # Add a red line to show which  is the current cutoff.
      points(x = input$mInput, y = 0, pch = "|", col="red", cex = 1.5)
    }, bg = "transparent")
    
    #' Generate the slider for the manual cutoff
    output$numericCutoff <- renderUI({
      sliderInput(inputId = "mInput",label = NULL, min = min(mRNAsurv()), max = max(mRNAsurv()), 
                  value = median(mRNAsurv()), step = 0.05, round = -2)
    })
    
    #' Create a Kaplan Meier plot with cutoff based on quantiles or manual selection
    output$survPlot <- renderPlot({     
      survNeed ()
      validate(need(input$histologySurv %in% c("All", histo()),""))   
      # Use 'try' to suppress a message throwed the first time manual cutoff is selected
      if(allSubSurv()) {
        try({
          par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2.2,.95,0))
          for (i in c("Classical","Mesenchymal","Neural","Proneural")) {
            survivalPlot (exprs(), input$gene, group = "GBM", cutoff = input$cutoff, numeric = input$mInput,
                          subtype = i, gcimp = gcimpSurv(), primary = primarySurv())}
        }, silent = TRUE)} else {
          try(survivalPlot (exprs(), input$gene, group = input$histologySurv, cutoff = input$cutoff, numeric = input$mInput,
                            subtype = input$subtypeSurv, gcimp = gcimpSurv(), primary = primarySurv()), silent = TRUE)}
    }, height = function(){if(!allSubSurv()) {400} else {600}}, width = function(){if(!allSubSurv()) {500} else {800}})
    
    
    #' Download the survPlot
    output$downloadsurvPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_survPlot.pdf")
      },      
      content = function(file) {
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(file, width = downloadPlotWidth(), height = downloadPlotHeight())
        if(allSubSurv()) {
          par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2.2,.95,0))
          for (i in c("Classical","Mesenchymal","Neural","Proneural")) {
            survivalPlot (exprs(), input$gene, group = "GBM", cutoff = input$cutoff, numeric = input$mInput,
                          subtype = i, gcimp = gcimpSurv(), primary = primarySurv())}
        } else {
          survivalPlot (exprs(), input$gene, group = input$histologySurv, cutoff = input$cutoff, numeric = input$mInput,
                        subtype = input$subtypeSurv, gcimp = gcimpSurv(), primary = primarySurv())}
        dev.off()
      }
    ) 
    
    #' Download the kmPlot
    output$downloadkmPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_kmPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(file, width = downloadPlotWidth(), height = downloadPlotHeight())
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
    output$corrPlot <- renderPlot({    
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        need(input$gene2 != "", "Please enter Gene 2")%then%
          need(input$gene2 %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the choosen histology is not available.
        need(input$histologyCorr %in% c("All",histo()),"")
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
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_corrPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(file, width = downloadPlotWidth(), height = downloadPlotHeight())
        myCorggPlot (exprs(), input$gene, input$gene2, input$histologyCorr, input$subtype, 
                     colorByInput(), separateByInput())
        dev.off()
      }
    )
    
    #' Multiple genes correlation
    pairsData <- reactive({
      validate(
        # Need two or more genes
        need(length(input$genelist) > 1, "Please enter 2 or more genes in the panel on the left")%then%
          need(input$genelist %in% names(exprs()),"Gene not available for this dataset")%then%         
          # Trying to avoid an error when switching datasets in case the choosen histology is not available.
          need(input$histologyCorr %in% c("All",histo()),"")
      )
      df <- exprs()
      if (input$histologyCorr != "All") {
        df <- subset (df, Histology == input$histologyCorr)
      } 
      if (input$histologyCorr == "GBM" & input$subtype != "All") {
        df <- subset (df, Subtype == input$subtype)
      }
      df <- df[ ,input$genelist] 
    })
    
    #' Generate the pairs plot
    output$pairsPlot <- renderPlot({
      #       myPairsPlot(pairsData())
      #       theme_set(theme_bw())
      ggpairs(pairsData(),lower=list(continuous="smooth", params=list(alpha=0.5)))  
    })
    
    #     #' Generate an HTML table view of the correlation table 
    #     output$pairsData <- renderTable({
    #       if (length(input$genelist) < 2)
    #         return()
    #       pairs.table <- rcorr(as.matrix(pairsData()))
    #       pairs.table
    #     })
    
    #' Download the pairs plot
    output$downloadpairsPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$dataset, "_pairsPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(downloadPlotFileType())
        plotFunction(file, width = downloadPlotWidth(), height = downloadPlotHeight())
        print(ggpairs(pairsData(),lower=list(continuous="smooth", params=c(alpha=0.5)))) 
        dev.off()
      }
    )
    
    #' Reactive function for an HTML table view of the data
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
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, ".csv")
      },
      content = function(file) {
        write.csv(dataTable(), file)
      }
    )
    
    #' Generate a graphic summary of the dataset, using ggvis
    output$piePlots <- renderUI({
      data <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status")]
      data <- data[,colSums(is.na(data)) < nrow(data)] # Removing unavailable (all NA) groups
      plot_output_list <- lapply(names(data), function(i) {
        plotname <- paste("plot", i, sep="")
        htmlOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items to display properly.
      do.call(tagList, plot_output_list)
    })
    
    observe ({                                                               
      data <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status")]
      data <- data[ ,colSums(is.na(data)) < nrow(data)]
      # Call renderChart for each one. 
      for (i in names(data)) {                                                    
        local({
          my_i <-i
          plotname <- paste("plot", my_i, sep="")
          output[[plotname]] <- renderGvis({
            plotData <- data.frame(table(data[, my_i]))
            pie <-gvisPieChart(labelvar = "Var1", numvar = "Freq", data = plotData, 
                               options = list(width=400, height=300, pieSliceText='label', chartArea.left = 1), chartid= plotname) # title = my_i,
            return(pie)
          })
        })
      }
    })
    
    #' Generate survival groups stratified by Histology, etc.
    output$survPlots <- renderUI({
      validate(
        need(!input$dataset %in% noSurvDataset,
             "Sorry, no survival data are available for this dataset")
      )
      df <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status","survival", "status")]
      df <- df[,colSums(is.na(df)) < nrow(df)] # Removing unavailable (all NA) groups
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      plot_output_list <- lapply(groups, function(i) {
        plot_surv_name <- paste("plotSurv", i, sep = "")
        plotOutput(plot_surv_name, height = 300, width = 400)
      })
      do.call(tagList, plot_output_list)
    })  
    
    observe({   
      df <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status","survival", "status")]
      df <- df[,colSums(is.na(df)) < nrow(df)] 
      groups <- names(df)[!names(df) %in% c("Sample","status","survival")]
      for (i in groups) {                                                    
        local({
          my_Survi <- i
          plot_surv_name <- paste("plotSurv", my_Survi, sep="")
          output[[plot_surv_name]] <- renderPlot({
            df1 <- na.omit(data.frame(status = df[ ,"status"], time = df[ ,"survival"], strata = df[ ,my_Survi]))
            df1$strata <- droplevels(df1$strata)
            expr.surv <- survfit(Surv(time, status == 1) ~ strata, data = df1, conf.type = "none")
            plot(expr.surv, xlab = "Survival time (Months)", ylab = "% Surviving", 
                 yscale = 100,  col = 1:length(df1$strata), mark.time = FALSE, main = paste(my_Survi))
            legend("topright", legend = levels(df1$strata), col = 1:length((levels(df1$strata))), lty = 1)
          })
        })
      }
    })
    
    #' Generate reports
    output$reportPlots <- renderUI({
      data <- data()
      groups <- c(plotList[[input$dataset]], plotUserSelection())
      plot_output_list <- lapply(groups, function(i) {
        plot_report <- paste("plotReport", i, sep = "")
        plotOutput(plot_report, height = 300)
      })
      do.call(tagList, plot_output_list)
    })  
      
    observe({
      data <- data()
      groups <- c(plotList[[input$dataset]], plotUserSelection())
      for (i in groups) {
        local({
          my_i <- i
          plot_report <- paste("plotReport", my_i, sep = "")
          output[[plot_report]] <- renderPlot({
            data <- subset(data,!is.na(data[,my_i]))
            p <- ggplot(data, mapping=aes_string(x=my_i, y = "mRNA")) + geom_boxplot(outlier.size = 0) +  
              geom_jitter(position = position_jitter(width = .2), size = 2, alpha = 0.5) + theme_bw() +
              ylab("mRNA expression (log2)") + xlab(paste0("\n",my_i)) + theme(axis.title.y=element_text(vjust=1)) 
            print(p)
          })
        })
      }
    })
    
    #' Generate reports table
    output$reportTables <- renderUI({
      data <- data()
      groups <- c(plotList[[input$dataset]], plotUserSelection())
      plot_output_list <- lapply(groups, function(i) {
        plot_report_summary <- paste("plotTable", i, sep = "")
        box(height = 280, title = paste0(i), width = NULL, solidHeader = TRUE, color = "primary",
        tableOutput(plot_report_summary)
        )
      })
      do.call(tagList, plot_output_list)
    }) 
    
    observe({
      groups <- c(plotList[[input$dataset]], plotUserSelection())
      for (i in groups) {
        local({
          my_i <- i
          plot_report_summary <- paste("plotTable", my_i, sep = "")
          output[[plot_report_summary]] <- renderTable({ 
            data <- data()
            validate(need(my_i %in% names(data),""))
            data <- subset(data,!is.na(data[,my_i]))
            stat <- substitute(data %>%
                                 group_by(x) %>%
                                 summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data())[1],2), "%)" ), # prop.table
                                           median = median (mRNA, na.rm=T), mad = mad(mRNA, na.rm=T),mean = mean(mRNA, na.rm=T), 
                                           sd = sd(mRNA, na.rm=T)),
                                 list(x = as.name(my_i)))
            stat <- data.frame(eval(stat))
            row.names(stat) <- stat[,my_i]
            stat <- stat[,-1]
            stat      
          })
        })
      }
    })
    
    #' Reactivity required to display download button after file upload
    output$finishedUploading <- reactive({
      if (is.null(input$upFile))
      { 0 } else { 1 }
    })
    outputOptions(output, 'finishedUploading', suspendWhenHidden=FALSE)
    
    #' Reactive function to generate SVM subtype call to pass to data table and download handler
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
        train <- train[train$Sample %in% gbm.core.samples,]
      } else if (input$tumorType == "lgg") {
        train <- lgg.tcga[["expr"]]
        train <- train[train$Sample %in% lgg.core.samples,]
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
      Training <- train$Subtype
      svm <- ksvm(t(df.train), Training, cross=10, kernel="vanilladot", family="multinomial", prob.model=TRUE, scale=FALSE)  
      svm.subtype.call <- as.matrix(predict(svm, t(df.learn)))
      if (input$tumorType == "gbm") {
        svm.subtype.call <- factor(svm.subtype.call,levels = c("Classical", "Mesenchymal", "Neural", "Proneural"))
      }
      prob <- as.matrix(predict(svm, t(df.learn), type="probabilities"))
      svm.call <- data.frame(Sample = rownames(upData), svm.subtype.call, prob)
      svm.call
    })
    
    #' Rerndering the subtype call as a data table
    output$svm <- renderDataTable({ 
      if (is.null(input$upFile) || input$goSvm == 0)
        return(NULL)
      input$goSvm
      isolate({
        svm <- svm.call()
      })
    })
    
    #' Download the subtype call
    output$downloadSvm <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", "SVM_Subtype_Call.csv")
      },
      content = function(file) {
        write.csv(svm.call(), file)
      }
    )
    
    #' Reactive function to generate k-nearest neighbour subtype call to pass to data table and download handler
    knn.call <- reactive ({
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      if (input$tumorType == "gbm") {
        train <- gbm.tcga[["expr"]]
        train <- train[train$Sample %in% gbm.core.samples,]
      } else if (input$tumorType == "lgg") {
        train <- lgg.tcga[["expr"]]
        train <- train[train$Sample %in% lgg.core.samples,]
      }
      row.names(train) <- train[,"Sample"]
      train.exp <- as.matrix(train[,-c(1:8)])
      row.names(upData) <- upData[,"Sample"]
      learn.exp <-as.matrix(upData [,-1])
      # Common genes of the two datasets
      genes <- intersect(colnames(train.exp), colnames(learn.exp))
      pred <- knn(train = train.exp[,genes], test = learn.exp[,genes], cl =  train$Subtype, k = 4, prob=TRUE)
      kn <- data.frame(Sample = rownames(upData), knn.subtype.call = pred, prob = attr(pred,"prob"))
      kn
    })
    
    #' Rerndering the knn subtype call as a data table
    output$knn <- renderDataTable({ 
      if (is.null(input$upFile) || input$goKnn == 0)
        return(NULL)
      input$goknn
      isolate({
        knn <- knn.call()
      })
    })
    
    #' Download the knn subtype call
    output$downloadKnn <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", "KNN_Subtype_Call.csv")
      },
      content = function(file) {
        write.csv(knn.call(), file)
      }
    )  
    #' Reactive function to generate ssGSEA call to pass to data table and download handler
    gsva.call <- reactive ({
      validate(
        need(input$tumorType == "gbm","ssGSEA analysis currently available for GBM samples only"))
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      if (input$tumorType == "gbm") {
        gene_list <- gbm.subtype.list
      } else if (input$tumorType == "lgg") {
        return(NULL)
      }
      rownames(upData) <- upData$Sample
      exprs <- data.frame(t(upData[,-1]))
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq=FALSE,
                           min.sz=0, max.sz=10000, verbose=FALSE)
      subtype_scores <- t(gsva_results)
      subtype_final <- data.frame(Sample = rownames(upData), gsea.subtype.call = names(gbm.subtype.list)[apply(subtype_scores,1,which.max)], 
                                  subtype_scores) 
      subtype_final
    })
    
    #' Rerndering the subtype call as a data table
    output$gsva <- renderDataTable({ 
      if (is.null(input$upFile) || input$goGsva == 0)
        return(NULL)
      input$goGsva
      isolate({
        gsva <- gsva.call()
      })
    })
    
    #' Download the subtype call
    output$downloadGsva <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", "GSVA_Subtype_call.csv")
      },
      content = function(file) {
        write.csv(gvsa.call(), file)
      }
    )
    
    #' Reactive function to generate the 3 subtype calls to pass to data table and download handler
    sub3.call <- reactive ({
      validate(
        need(input$tumorType == "gbm","ssGSEA analysis currently available for GBM samples only"))
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
      sub3 <- data_frame(Sample = upData$Sample, svm.call = svm.call()[,"svm.subtype.call"], 
                         knn.call = knn.call()[,"knn.subtype.call"], gsea.call = gsva.call()[,"gsea.subtype.call"],
                         equal.call = ifelse(svm.call==knn.call & knn.call==gsea.call, "TRUE", "FALSE"))
      sub3$majority.call <- apply(sub3[,2:4],1,maj) # compare calls
      sub3 <- as.data.frame(sub3) # the dataframe generated with dplyr gives issue when sorting the html table
    })
    
    #' Rerndering the subtype call as a data table
    output$sub3 <- renderDataTable({ 
      if (is.null(input$upFile) || input$goSub3 == 0)
        return(NULL)
      input$goSub3
      isolate({
        sub3 <- sub3.call()
      })
    })
    
    #' Download the subtype call
    output$downloadSub3 <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", "3way_Subtype_call.csv")
      },
      content = function(file) {
        write.csv(sub3.call(), file)
      }
    )
    
    #' Correlation method
    corrMethod <- reactive({
      switch(input$corrMethod,
             "Spearman" = "spearman",
             "Pearson" = "pearson")
    })
    
    #' Generate the correlation table ##  corFast works locally but not on shinyapps.io
    corr <- reactive ({
      corr <- getCorr(datasetInputCor()[["expr"]], input$geneCor, input$histologyCorrTable, corrMethod())
      corr  <- merge(genes, corr, by="Gene")
#       corr <- corr [-1,]
      corr <- corr[order(-abs(corr$r)), ]
    })
    
    #' Generate a reactive element of the the correlation data 
    corrData <- reactive({  
      corr.table <- suppressWarnings(corr())  # suppressWarnings  is used to prevent the warning messages in the LGG dataset  
      corr.table <- subset(corr.table, adj.p.value <= as.numeric(input$sign))
      if (input$cor == "Positive"){
        corr.table <- subset(corr.table, r > 0)
        corr.table <- corr.table[order(-corr.table$r), ]
      }
      if (input$cor == "Negative"){
        corr.table <- subset(corr.table, r < 0)
        corr.table <- corr.table[order(corr.table$r), ]
      } 
      if (input$cor == "All"){
        corr.table <- subset(corr.table, r <= input$range[1] | r >= input$range[2])
      }
      corr.table
    })
    
    #' Generate an HTML table view of the correlation table 
    output$corrData <- renderDataTable({
      validate(
        need(input$geneCor != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$geneCor %in% names(datasetInputCor()[["expr"]]),"Gene not available for this dataset")
      )   
      corrData()
    }, options = list(orderClasses = TRUE), callback = "function(table) {
      table.on('click.dt', 'tr', function() {
            table.$('tr.selected').removeClass('selected');
            $(this).toggleClass('selected');            
        Shiny.onInputChange('rows',
                            table.rows('.selected').data()[0][0]);
      });
    }")
    
    #' Generate a reactive value for the input$rows that set to NULL when the dataset change
    v <- reactiveValues(rows = NULL)
    observeEvent(input$rows, {
      v$rows <- input$rows
    })
    observeEvent(datasetInputCor(), {
      v$rows <- NULL
    })
    observeEvent(input$histologyCorrTable, {
      v$rows <- NULL
    })
    
    #' Generate the correlation plot
    output$corrDataPlot <- renderPlot({
      v$rows
      validate(
        need(input$geneCor != "","")%then%
          need(v$rows!= "","Click on a row to see the corresponding correlation plot.")
      )
      df <- datasetInputCor()[["expr"]]
      if (input$histologyCorrTable != "All") {
        df <- subset (df, Histology == input$histologyCorrTable)
      } else {
        df <- df
      }
      aes_scatter <- aes_string(x = input$geneCor, y = v$rows)
      ggplot(df,mapping = aes_scatter) + theme(legend.position=c(1,1),legend.justification=c(1,1)) +
        geom_point(alpha=.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1)
      
    })
    
    
    #' Download the correlation table 
    output$downloadCorrData <- downloadHandler(
      filename = function() {
        paste(Sys.Date(), input$geneCor, input$datasetCor, input$histologyCorrTable, 
              "corrData.csv", sep="_")
      },
      content = function(file) {
        write.csv(corrData(),file)
      }
    )
    
  })
