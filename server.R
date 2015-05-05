# server.R for Gliovis
shinyServer(
  function(input, output, session) {
    
    options(shiny.maxRequestSize=30*1024^2)    
    
    observe({
      if (input$nav == "quit") stopApp()
    })
    
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
    
    #' RPPA data
    rppas <- reactive({
      datasetInput()[["rppa"]]
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
      # Exclude the pre-defeined plots and numeric variabe
      drop <- c("Sample","Histology","Grade","Recurrence","Subtype", "CIMP_status", "survival",
                "status", "Age", "ID","Patient_ID","Sample_ID", "Matching.sample", "Therapy_Class","title") 
      data <- pDatas()[,!names(pDatas()) %in% drop]
      n <- names(data)
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
    
    plotType <- reactive({
      if (input$plotType == "Pre-defined"){
        plotSelected()
      } else if (input$plotType == "User-defined") {
        plotUserSelected()
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
      data <- cbind(mRNA, exprs()[,2:6]) # To combine with pData
      data <- cbind(data, pDatas()[,!names(pDatas()) %in% names(data)]) # To combine with more pData for the report
      if (input$dataset == "TCGA GBM" | input$dataset == "TCGA Lgg") {
        if(input$gene %in% names(cnas())){
          Copy_number <- cnas()[ ,input$gene]}
        else {
          Copy_number <- rep(NA, length(cnas())) # Some genes don't have copy numbers data
        }
        Copy_number <- factor(Copy_number, levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp")) 
        data <- cbind(Copy_number,data)
      }
      data 
    })
       
    #' Data for the box plot
    plotData <- reactive({
      validate(
        need(!all(is.na(data()[ ,plotType()])),"Sorry,no gene data available for this group")
      )
      data <- data()
      if (input$primary & any(!is.na(data$Recurrence))) {
        data <- filter(data, Recurrence == "Primary")
      }
      data <- subset(data, !is.na(data[ ,plotType()]))
    })
    
    #' Reactive function to generate the box plots
    boxplot <- reactive ({
      data <- plotData()
      if (input$scale) {
        ylab <- "Normalized mRNA expression"
      } else {
        ylab <- "mRNA expression (log2)"
      }
      if (input$colBox) {
        box <- geom_boxplot(aes_string(fill = plotType()), outlier.size = 0) # It works but not the right way to approach this issue
      } else {
        box <- geom_boxplot(outlier.size = 0)
      }
      if (input$colStrip) {
        col <- aes_string(color = input$colorP)
        strip <- geom_jitter(position = position_jitter(width = .2), col, size = 2, alpha = 0.75)
      } else {
        strip <- geom_jitter(position = position_jitter(width = .2), size = 2, alpha = 0.5)
      }
      p <- ggplot(data, mapping=aes_string(x=plotType(), y = "mRNA")) + ylab(ylab) + xlab(paste0("\n",plotType())) +
        theme(axis.title.y=element_text(vjust=1)) + theme_gray(base_size = 14)
      p <- p + box + strip
      if (input$bw) {
        p <- p + theme_bw (base_size = 14) 
      }
      if (input$tukeyPlot) {
        mRNA <- data[,"mRNA"]
        group <- data[ ,plotType()]
        tukey <- data.frame(TukeyHSD(aov(mRNA ~ group))[[1]])
        t <- tukey %>%
          mutate(Significance = as.factor(starmaker(p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns"))),
                 comparison = row.names(.)) %>%
          ggplot(aes(reorder(comparison, diff), diff, ymin = lwr, ymax= upr, colour = Significance)) +
          geom_point() + geom_errorbar(width = 0.25) + ylab("Differences in mean levels") + xlab("") + 
          geom_hline(xintercept = 0, colour="darkgray", linetype = "longdash") + coord_flip()
        if (input$bw) {
          t <- t + theme_bw ()
        }
        grid.arrange(p, t, ncol=2, widths = c(3,2))
      } else {
        p
      }
    })
  
    #' Generate radiobuttons for the various categories in the pData
    output$colorPoints <- renderUI({
      validate(
        need(input$gene != "","")
      )
      data <- rmNA(plotData())
      colnames <- names(data)[!names(data) %in% c("mRNA","Sample","status","survival")]
      # Create the radiobuttons for the different pData categories
      radioButtons("colorP", "Color by:", choices  = colnames, selected = plotType())
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({
      # To avoid an error when switching datasets in case the colStrip is not available.
      if(input$colStrip){
        data <- rmNA(plotData())
        colnames <- names(data)[!names(data) %in% c("mRNA","Sample","status","survival")]
        validate(need(input$colorP %in% colnames,""))
      }
      print(boxplot())
    })
    
    #' Summary statistic
    output$summary <- renderTable({    
      data <- plotData()
      data$group <- data[ ,plotType()]
      stat <- data.frame(data %>%
                           group_by(group) %>%
                           summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data)[1],2), "%)" ), # prop.table
                                     median = median (mRNA, na.rm=T), mad = mad(mRNA, na.rm=T),mean = mean(mRNA, na.rm=T), sd = sd(mRNA, na.rm=T)))
      row.names(stat) <- stat$group
      tot <- data %>%
        summarise(Sample_count = n(),median = median (mRNA, na.rm=T), 
                  mad = mad(mRNA, na.rm=T),mean = mean(mRNA, na.rm=T), sd = sd(mRNA, na.rm=T))
      stat <- stat[,-1]
      stat <- rbind(stat,TOTAL = tot)
      stat 
    }, align='rrrrrr')

    #' Tukey post-hoc test
    output$tukeyTest <- renderTable({    
      data <- plotData()
      mRNA <- data[,"mRNA"]
      group <- data[ ,plotType()]
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff), ]
      tukey
    })
    
    #' Pairwise t test
    output$pairwiseTtest <- renderTable({     
      data <- plotData()
      mRNA <- data[,"mRNA"]
      group <- data[ ,plotType()]
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
        print(boxplot())  
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
    
    #' Extract the survival data.
    survData <- reactive({    
      df <- data()
      df <- filter(df, !is.na(df$status))
      if (input$histologySurv != "All"){
        df <- filter(df, Histology == input$histologySurv)
      }
      if (input$histologySurv == "GBM" & input$subtypeSurv != "All") {
        df <- filter (df, Subtype == input$subtypeSurv)
      }
      # exclude G-CIMP is selected
      if (gcimpSurv()){
        df <- filter (df, CIMP_status != "G-CIMP")
      }
      # select primary sample
      if (primarySurv() & any(!is.na(df$Recurrence))) {
        df <- filter (df, Recurrence == "Primary")
      }
      df
    })  

    #' Subset to GBM samples for the interactive HR plot.
    survGBM <- reactive({
      df <- filter (survData(), Histology == "GBM")
    })
    
    #' Extract the GBM expression values for the interactive HR plot.
    geneExp <- reactive({
      geneExp <- survGBM()[ ,"mRNA"]
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
      HR <- getHR(survGBM())
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
      validate(need(!input$dataset %in% c("TCGA Lgg","Gorovets"), "Interactive HR plot currently available only for GBM samples") %then%
                 need(histoSurvSelected() == "GBM","Please select GBM samples in the 'Histology' dropdown menu") %then%
                   need(!input$dataset %in% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      survNeed()
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
      surv <- with(survGBM(), Surv(survival, status == 1))
      return(surv ~ expressionGrp)
    })
    
    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      validate(need(histoSurvSelected() == "GBM", ""))
      survNeed()
      cutoff <- getCutoff()
      surv <- survivalFml()
      kmPlot(cutoff, surv)
    })
    
    #' Create a slider for the manual cutoff of the Kaplan Meier plot
    mRNAsurv <- reactive({
      survNeed()
      validate(need(input$histologySurv %in% c("All", histo()),""))      
      mRNA <- survData()[ ,"mRNA"]
      mRNA.values <- round(mRNA[!is.na(mRNA)],2)
      # Generate a vector of continuos values, excluding the first an last value
      mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) 
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxmRNA <- renderPlot({    
      validate(need(input$mInput, ""))      
      mRNA <- round(mRNAsurv(),2)
      q <- quantile(mRNA)
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|")
      # Add a red line to show which  is the current cutoff.
      points(x = input$mInput, y = 0, pch = "|", col="red", cex = 2.5)
      points(x = q[2:4], y = rep(0,3), pch = "|", col="blue", cex = 2)
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
            survivalPlot (survData(), input$gene, group = "GBM", subtype = i,
                          cutoff = input$cutoff, numeric = input$mInput)
          }
        }, silent = TRUE)} else {
          try(survivalPlot (survData(), input$gene, group = input$histologySurv, subtype = input$subtypeSurv,
                            cutoff = input$cutoff, numeric = input$mInput), silent = TRUE)
        }
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
            survivalPlot (survData(), input$gene, group = "GBM", subtype = i,
                          cutoff = input$cutoff, numeric = input$mInput)
          }
        } else {
          survivalPlot (survData(), input$gene, group = input$histologySurv, subtype = input$subtypeSurv,
                        cutoff = input$cutoff, numeric = input$mInput)
        }
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
        df <- filter(df, Histology == input$histologyCorr)
      } 
      if (input$histologyCorr == "GBM" & input$subtype != "All") {
        df <- filter(df, Subtype == input$subtype)
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
     
    #' RPPA data analysis
    rppaRNA <- reactive({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA Lgg"), "")%then%
          need(input$gene != "", "")%then%
          need(input$gene %in% names(exprs()),"")
      )
      samples <- intersect(row.names(rppas()),exprs()[,"Sample"])
      mRNA <- round(exprs()[samples,input$gene],2)
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxRppaRNA <- renderPlot({    
      validate(need(input$rppaCut, ""))      
      mRNA <- rppaRNA()
      q <- quantile(mRNA)
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|")
      # Add a red line to show which  is the current cutoff.
      points(x = input$rppaCut, y = 0, pch = "|", col="red", cex = 2)
      points(x = q[2:4], y = rep(0,3), pch = "|", col="blue", cex = 2)
      #       abline(v= q[2:4], col="blue")
    }, bg = "transparent")
    
    output$rppaCutoff <- renderUI({
      sliderInput(inputId = "rppaCut",label = "mRNA cutoff", min = min(rppaRNA()), max = max(rppaRNA()), 
                  value = median(rppaRNA()), step = 0.05, round = -2)
    })
       
    output$rppaTable <- renderDataTable({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA Lgg"), "RPPA data available only for TCGA datasets")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")%then%
          need(input$rppaCut != "", "")%then%
          need(input$rppaCut > min(rppaRNA()) & input$rppaCut < max(rppaRNA()), "")
      )
      rppa <- rppas()
      samples <- intersect(row.names(rppa),exprs()[,"Sample"])
      mRNA <- round(exprs()[samples,input$gene],2)
      strat <- ifelse(mRNA >= input$rppaCut, c("high"),c("low"))
      strat <- factor(strat,levels = c("low", "high"))
      d <- NULL
      results <- for (i in names(rppa)) {
        prot <- rppa[, i]
        t <- t.test(prot~strat)
        est <- t$estimate
        p <-  t$p.value
        x <- c(est,p)
        names(x) <- c("Avg_low","Avg_high","p.value")
        d= rbind(d,t(data.frame(x)))
      }
      row.names(d) <- names(rppa)
      d <- data.frame(Protein = row.names(d), round(d,5))
      # d$adj.p.value <- p.adjust(d$p, method = "bonferroni")
      d <- d[order(d$p.value),]
      datatable(d, rownames = FALSE, options = list(lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full"),
                callback = JS("
                table.on('click.dt', 'tr', function() {
                table.$('tr.selected').removeClass('selected');
                $(this).toggleClass('selected');            
                Shiny.onInputChange('rppa.rows',
                table.rows('.selected').data()[0][0]);
                })")
      )
      })
    
    #' Generate a reactive value for the input$rows that set to NULL when the dataset change
    rp <- reactiveValues(rppa.rows = NULL)
    observeEvent(input$rppa.rows, {
      rp$rppa.rows <- input$rppa.rows
    })
    observeEvent(datasetInput(), {
      rp$rppa.rows <- NULL
    })
    observeEvent(input$histology, {
      rp$rppa.rows <- NULL
    })
    
    observeEvent(input$gene, {
      rp$rppa.rows <- NULL
    })
    
    #' Generate the RPPA box plot
    output$rppaPlot <- renderPlot({
      rp$rppa.rows
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA Lgg"), "")%then%
          need(input$gene != "","")%then%
          need(rp$rppa.rows!= "","Click on a row to see the corresponding plots.")
      )
      mRNA <- rppaRNA()
      strat <- ifelse(mRNA >= input$rppaCut, c("high"),c("low"))
      strat <- factor(strat,levels = c("low", "high"))
      data <- data.frame(mRNA, strat, rppa = rppas()[,rp$rppa.rows])
      p1 <- ggplot(data, aes(x=strat, y = rppa)) + geom_boxplot(outlier.size = 0) + 
        geom_jitter(aes(colour = strat), position = position_jitter(width = .2), size = 2, alpha = 0.5) + 
        xlab(paste(input$gene, "mRNA")) + ylab(paste(rp$rppa.rows,"RPPA score")) + 
        guides(colour=FALSE) + theme_bw()
      p2 <- ggplot(data, aes(x=mRNA, y = rppa)) + geom_point(aes(colour = strat), alpha=.5) +
        geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_bw() +
        xlab(paste(input$gene, "mRNA (log2)")) + ylab(paste(rp$rppa.rows,"RPPA score")) + theme(legend.position = "none")
      grid.arrange(p1, p2, ncol=1)
    },height = 700)
    
    #' Generate reports
    output$reportPlots <- renderUI({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")
      )
      groups <- c(plotList[[input$dataset]], plotUserSelection())
      plot_output_list <- lapply(groups, function(i) {
        plot_report <- paste("plotReport", i, sep = "")
        box(height = 300, title = paste0(i), width = NULL, solidHeader = TRUE, status = "primary",
            plotOutput(plot_report, height = 245)
        )
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
            validate(
              need(!all(is.na(data[ ,my_i])),"Sorry,no gene data available for this group")
            )
            data <- filter(data,!is.na(data[,my_i]))
            p <- ggplot(data, mapping=aes_string(x=my_i, y = "mRNA")) + geom_boxplot(outlier.size = 0) +  
              geom_jitter(position = position_jitter(width = .2), size = 2, alpha = 0.5) + 
              ylab("mRNA expression (log2)") + theme(axis.title.x = element_blank()) + theme(axis.title.y=element_text(vjust=1)) 
            stat <- substitute(data %>%
                                 group_by(x) %>%
                                 summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data())[1],2), "%)" ), # prop.table
                                           median = round(median (mRNA, na.rm=T),2), mad = round(mad(mRNA, na.rm=T),2),
                                           mean = round(mean(mRNA, na.rm=T),2), sd = round(sd(mRNA, na.rm=T),2)),
                               list(x = as.name(my_i)))
            stat <- data.frame(eval(stat))
            row.names(stat) <- stat[,my_i]
            tot <- eval(substitute(data %>%
                                     summarise(Sample_count = n(), median = round(median (mRNA, na.rm=T),2), mad = round(mad(mRNA, na.rm=T),2),
                                               mean = round(mean(mRNA, na.rm=T),2), sd = round(sd(mRNA, na.rm=T),2)),
                                   list(x = as.name(my_i))))
            stat <- stat[,-1]
            stat <- rbind(stat,TOTAL = tot)
            t <- tableGrob(stat, gp = gpar(fontsize=14),row.just = "right", core.just = "right")
            grid.arrange(p, t, ncol = 2, just = c("center", "top")) # `just` it's not working
          })
        })
      }
    })
    
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
      data <- rmNA(data)
      data
    })
    
    #' Generate an HTML table view of the data
    output$table <- renderDataTable({
      # If a gene is not specified show the pData only
      if (input$gene == "") {
        data <- rmNA(pDatas())
      } else {
        data <- dataTable()
      }
      DT::datatable(data, rownames = FALSE, extensions = c("FixedColumns", "TableTools"),
                    options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE,
                                   lengthMenu = c(10, 30, 50), pageLength = 10, dom = 'T<"clear">lfrtip',
                                   tableTools = list(aButtons = c("copy","csv","xls","print"),
                                                     sSwfPath = copySWF(dest = "www"))))
    })
    
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
      df <- rmNA(df) 
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
    
    #' Correlation method
    corrMethod <- reactive({
      switch(input$corrMethod,
             "Spearman" = "spearman",
             "Pearson" = "pearson")
    })
    
    #' Generate the correlation table 
    corr <- reactive ({
      corr <- getCorr(datasetInputCor()[["expr"]], input$geneCor, input$histologyCorrTable, corrMethod())
      corr  <- merge(genes, corr, by="Gene")
      corr <- arrange(corr, desc(abs(r)))
    })
    
    #' Generate a reactive element of the the correlation data 
    corrData <- reactive({  
      corr.table <- suppressWarnings(corr())  # suppressWarnings  is used to prevent the warning messages in the LGG dataset  
      corr.table <- filter(corr.table, adj.p.value <= as.numeric(input$sign))
      if (input$cor == "Positive"){
        corr.table <- filter(corr.table, r > 0)
        corr.table <- arrange(corr.table, desc(r))
      }
      if (input$cor == "Negative"){
        corr.table <- filter(corr.table, r < 0)
        corr.table <- arrange(corr.table,r)
      } 
      if (input$cor == "All"){
        corr.table <- filter(corr.table, r <= input$range[1] | r >= input$range[2])
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
      DT::datatable(corrData(), rownames = FALSE, extensions = "TableTools",
                    options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full",
                                   dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                               sSwfPath = copySWF(dest = "www"))),
                    callback = JS("
                      table.on('click.dt', 'tr', function() {
                        table.$('tr.selected').removeClass('selected');
                        $(this).toggleClass('selected');
                        Shiny.onInputChange('rows', table.rows('.selected').data()[0][0]);
                      })"
                    )
                    )
    })
    
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
    observeEvent(input$geneCor, {
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
        df <- filter (df, Histology == input$histologyCorrTable)
      } else {
        df <- df
      }
      aes_scatter <- aes_string(x = input$geneCor, y = v$rows)
      ggplot(df,mapping = aes_scatter) + theme(legend.position=c(1,1),legend.justification=c(1,1)) +
        geom_point(alpha=.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_bw() 
    })
    
    #' Reactivity required to display download button after file upload
    output$finishedUploading <- reactive({
      if (is.null(input$upFile))
      { 0 } else { 1 }
    })
    outputOptions(output, 'finishedUploading', suspendWhenHidden=FALSE)
    
    #' Reactive function to generate SVM subtype call to pass to data table and download handler
    #' IMPORTANT #http://stackoverflow.com/questions/15503027/why-are-probabilities-and-response-in-ksvm-in-r-not-consistent
    svm.call <- eventReactive (input$goSvm | input$goSub3,{
      # input$upFile will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,stringsAsFactors=FALSE)
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
      set.seed(1234)
      svm <- ksvm(t(df.train), Training, cross=10, kernel="vanilladot", family="multinomial", prob.model=TRUE, scale=FALSE)  
      svm.subtype.call <- as.matrix(predict(svm, t(df.learn)))
      if (input$tumorType == "gbm") {
        svm.subtype.call <- factor(svm.subtype.call,levels = c("Classical", "Mesenchymal", "Neural", "Proneural"))
      }
      prob <- as.matrix(predict(svm, t(df.learn), type="probabilities"))
      svm.call <- data.frame(Sample = rownames(upData), svm.subtype.call, round(prob,3))
      svm.call
    })
    
    #' Rerndering the subtype call as a data table
    output$svm <- renderDataTable({ 
      if (is.null(input$upFile) || input$goSvm == 0)
        return(NULL)
        DT::datatable(svm.call(), rownames = FALSE, extensions = "TableTools",
                      options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full",
                                     dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                                 sSwfPath = copySWF(dest = "www"))))
    })
    
    #' Reactive function to generate k-nearest neighbour subtype call to pass to data table and download handler
    knn.call <- eventReactive (input$goKnn | input$goSub3, {
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
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
      set.seed(1234)
      pred <- knn3Train(train = train.exp[,genes], test = learn.exp[,genes], cl =  train$Subtype, k = 15, prob=TRUE)
      kn <- data.frame(Sample = rownames(upData), knn.subtype.call = pred, prob = round(attr(pred,"prob"),2))
      names(kn)[3:6] <- c("Classical","Mesenchymal","Neural","Proneural")
      kn
    })
    
    #' Rerndering the knn subtype call as a data table
    output$knn <- renderDataTable({ 
      if (is.null(input$upFile) || input$goKnn == 0)
        return(NULL)
        DT::datatable(knn.call(), rownames = FALSE, extensions = "TableTools",
                      options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full",
                                     dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                                 sSwfPath = copySWF(dest = "www"))))
    })
    
    #' Reactive function to generate ssGSEA call to pass to data table and download handler
    gsva.call <- eventReactive(input$goGsva | input$goSub3, {
      validate(
        need(input$tumorType == "gbm","ssGSEA analysis currently available for GBM samples only"))
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,stringsAsFactors=FALSE)
      if (input$tumorType == "gbm") {
        gene_list <- gbm.subtype.list
      } else if (input$tumorType == "lgg") {
        return(NULL) # 
      }
      rownames(upData) <- upData$Sample
      exprs <- data.frame(t(upData[,-1]))
      set.seed(1234)
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq=FALSE, parallel.sz = 10,
                           min.sz=0, max.sz=10000, verbose=TRUE)
      subtype_scores <- round(t(gsva_results),3)
      subtype_final <- data.frame(Sample = rownames(upData), gsea.subtype.call = names(gene_list)[apply(subtype_scores,1,which.max)], 
                                  subtype_scores)
      subtype_final <- subtype_final[,c("Sample", "gsea.subtype.call", "Classical","Mesenchymal","Neural","Proneural")]
      subtype_final
    })
    
    #' Rerndering the subtype call as a data table
    output$gsva <- renderDataTable({ 
      if (is.null(input$upFile) || input$goGsva == 0)
        return(NULL)
        DT::datatable(gsva.call(), rownames = FALSE, extensions = "TableTools",
                      options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full",
                                     dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                                 sSwfPath = copySWF(dest = "www"))))
    })

    #' Reactive function to generate the 3 subtype calls to pass to data table and download handler
    sub3.call <- eventReactive (input$goSub3, {
      validate(
        need(input$tumorType == "gbm","'3-way' analysis currently available for GBM samples only"))
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,stringsAsFactors=FALSE)
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
        DT::datatable(sub3.call(), rownames = FALSE, extensions = "TableTools",
                      options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full", autoWidth = TRUE,
                                     dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                                 sSwfPath = copySWF(dest = "www"))))
    })
        
    #' Reactivity required to display download button after file upload
    output$finishedEstUploading <- reactive({
      if (is.null(input$upEstFile))
      { 0 } else { 1 }
    })
    outputOptions(output, 'finishedEstUploading', suspendWhenHidden=FALSE)
    
    #' Reactive function to generate Estimate call to pass to data table 
    est.call <- eventReactive (input$goEst,{
      inFile <- input$upEstFile
      upData <- read.csv(inFile$datapath, header=input$headerEst, sep=input$sepEst, quote=input$quoteEst, stringsAsFactors=FALSE)
      row.names(upData) <- upData[,"Sample"]
      ds <- data.frame(t(upData[,-1]))
      est <- myEstimateScore(ds, platform = input$platformEst)
      est <- cbind(Sample = upData$Sample, round(est,2))
    })
    
    #' Rerndering the subtype call as a data table
    output$estScore <- renderDataTable({ 
      if (is.null(input$upEstFile) || input$goEst == 0)
        return(NULL)
      DT::datatable(est.call(), rownames = FALSE, extensions = "TableTools", 
                    options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pageLength = 20, pagingType = "full", autoWidth = TRUE,
                                   dom = 'T<"clear">lfrtip', tableTools = list(aButtons = c("copy","csv","xls","print"), 
                                                                               sSwfPath = copySWF(dest = "www"))),
                    callback = JS("table.on('click.dt', 'tr', function() {
                                          table.$('tr.selected').removeClass('selected');
                                          $(this).toggleClass('selected');            
                                          Shiny.onInputChange('rowsEst',table.rows('.selected').data()[0][0]);
                                    })")
                    )
    })

    #' Generate the Purity plot
    output$purityPlot <- renderPlot({
      validate(
        need(est.call(),"")%then%
          need(input$platformEst == "affymetrix", "Sorry, the purity plot for Agilent and RNAseq have not yet been impemented") %then%
            need(input$rowsEst != "","Click on a row to see the corresponding purity plot.")       
      )
      plotPurity(est.call(), input$rowsEst, platform = input$platformEst)
    })
    
  })
