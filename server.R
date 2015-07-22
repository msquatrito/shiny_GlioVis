# server.R for Gliovis
shinyServer(
  function(input, output, session) {
    
    options(shiny.maxRequestSize=100*1024^2)    
    
    #' Return the requested dataset
    datasetInput <- reactive({
      switch(input$dataset, 
             "TCGA GBM" = gbm.tcga,
             "TCGA LGG" = lgg.tcga,
             "TCGA GBMLGG" = lgg_gbm.tcga,
             "Rembrandt" = rembrandt,
             "Gravendeel" = gravendeel,
             "Phillips" = phillips,
             "Murat" = murat,
             "Freije" = freije,
             "Lee Y" = leey,
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
             "Oh" = oh,
             "Ivy GAP" = ivy,
             "POLA Network" = pola,
             "Gleize" = gleize)
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
    # updateSelectizeInput(session, inputId = "geneCor", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "gene2", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "genelist", choices = gene_names, server = TRUE) 
    
    #' Required for the conditional panel 'corrMany' to work correctly
    observe({
      if(input$tab1 != 3)
        updateTabsetPanel(session, inputId = "tabCorr", selected = "corrTwo")
    })
    
    observe({
      if(input$tabTools != "DeconvoluteMe")
        updateCheckboxInput(session, inputId = "deconvPData",value = FALSE) 
    })
    
    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plot_Selected <- reactive ({
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
    plot_User_Selected <- reactive({
      if (input$plotTypeUserSel %in% plotUserSelection()){ 
        input$plotTypeUserSel
      } else {
        NULL
      }
    })
    
    plot_Type <- reactive({
      if (input$plotType == "Pre-defined"){
        plot_Selected()
      } else if (input$plotType == "User-defined") {
        plot_User_Selected()
      }
    })
    
    #' Change the plot type available for a specific dataset
    observe({
      updateSelectInput(session, inputId = "plotTypeSel", choices = plotList[[input$dataset]], selected = plot_Selected()) 
    }, priority = 10)
    
    observe({
      updateSelectInput(session, inputId = "plotTypeUserSel", choices = plotUserSelection(), selected = plot_User_Selected())
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
    histo_Surv_Selected <- reactive ({
      if (input$histologySurv %in% c("All", histo())){
        input$histologySurv
      } else {
        "All"
      }
    })
    
    #' When switching datasets for correlation plots, if the selected histo is not available it will choose "All"
    histo_Corr_Selected <- reactive ({
      if (input$histologyCorr %in% histo()){
        input$histologyCorr
      } else {
        return("All")
      }
    })
    
    observe({
      # This will change the value of input$histologySurv, based on histological group available for that dataset  
      updateSelectInput(session, inputId = "histologySurv", choices = c("All", histo()), selected = histo_Surv_Selected())
      # This will change the value of input$histologyCorr, based on histological group available for that dataset
      updateSelectInput(session, inputId = "histologyCorr", choices = c("All", histo()), selected = histo_Corr_Selected())
    })
    
    #' Generate a dataframe with the data to plot 
    data <- reactive({     
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the plotType is not available.
        need(plot_Selected() %in% plotList[[input$dataset]], FALSE),
        need(plot_User_Selected() %in% plotUserSelection(), FALSE)
      ) 
      mRNA <- exprs()[ ,input$gene]
      if (input$scale) {
        mRNA <- scale(mRNA)
      }
      data <- cbind(mRNA, exprs()[,2:6]) # To combine with pData
      data <- cbind(data, pDatas()[,!names(pDatas()) %in% names(data)]) # To combine with more pData for the report
      if (input$dataset %in% c("TCGA GBM", "TCGA LGG", "TCGA GBMLGG")) {
        if(input$gene %in% names(cnas())){
          Copy_number <- cnas()[ ,input$gene]}
        else {
          Copy_number <- rep(NA, length(cnas())) # Some genes don't have copy numbers data
        }
        Copy_number <- factor(Copy_number, levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
        Copy_number <- droplevels(Copy_number)
        data <- cbind(Copy_number,data)
      }
      data
    })
    
    #' Data for the box plot
    plot_Data <- reactive({
      validate(
        need(!all(is.na(data()[ ,plot_Type()])),"Sorry,no gene data available for this group")
      )
      data <- data()
      data <- subset(data, !is.na(data[ ,plot_Type()]))   
    })
    
    #' Filtered data for the box plot
    filter_plot_Data <- reactive({
      if (input$removeMe) {
        validate(need(input$removeGp %in% plot_Data()[ ,plot_Type()],message = FALSE)) #To silence an error thrown by the stat analysis
        data <- plot_Data()
        data <- subset(data, data[ ,plot_Type()] %in% input$removeGp)
      } else {
        data <- plot_Data()
      }
      data
    })
    
    #' Select input for groups to add/exclude
    output$removeGp <- renderUI({
      groups <- levels(plot_Data()[ ,plot_Type()])
      selectizeInput(inputId = "removeGp", label = "", choices = groups, selected = groups, multiple = TRUE, 
                     options = list(plugins = list('remove_button')))
    })
    
    #' Reactive function to generate the box plots
    box_Plot <- reactive ({
      data <- filter_plot_Data() 
      xlabel <- ifelse(input$labelsTitle, paste0("\n", input$myXlab), paste0("\n",plot_Type()))
      ylabel <- ifelse(input$scale, "Normalized mRNA expression\n", "mRNA expression (log2)\n")
      ylabel <- ifelse(input$labelsTitle, paste0("\n", input$myYlab), ylabel)
      theme <- theme(axis.text.x = element_text(size = input$axis_text_size), axis.text.y = element_text(size = input$axis_text_size),
                     axis.title.x = element_text(size = input$axis_title_size), axis.title.y = element_text(size = input$axis_title_size),
                     plot.margin = unit(c(0,0,0,0), "lines"))
      
      if (input$bw) {
        theme <- theme_bw () + theme
      }
      
      p <- ggplot(data, mapping=aes_string(x=plot_Type(), y = "mRNA")) + ylab(ylabel) + xlab(xlabel) + theme
     
       if (input$colBox) {
        p <- p + geom_boxplot(aes_string(fill = plot_Type()), outlier.size = 0)
      } else {
        p <- p + geom_boxplot(outlier.size = 0)
      }
      
      if (input$typePoint) {
        col <- input$colorP
        shape <- input$shapeP
        if(input$colorP == "None") {
          col <-  NULL
        }
        if(input$shapeP == "None") {
          shape <-  NULL
        }
        map <- aes_string(colour = col, shape = shape)
      } else {
        map <- NULL
      }
        
      p <- p + geom_jitter(position = position_jitter(width = .2), mapping = map, size = input$point_size, alpha = input$alpha)
      
      if (input$xaxisLabelAngle) {
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
      
      if (input$tukeyHSD && input$tukeyPlot) {
        t <- tukey() %>%
          mutate(comparison = row.names(.)) %>%
          ggplot(aes(reorder(comparison, diff), diff, ymin = lwr, ymax= upr, colour = Significance)) +
          geom_point() + geom_errorbar(width = 0.25) + ylab("\nDifferences in mean levels") + xlab("") + 
          geom_hline(xintercept = 0, colour="darkgray", linetype = "longdash") + coord_flip() + theme
        
        arrangeGrob(p, t, ncol=2, widths = c(3,2))
      } else {
        p
      }
    })
    
    #' Table with the data used for the plot
    output$filterDataTable <- renderDataTable({
      data <- filter_plot_Data()[,c("Sample", plot_Type(), "mRNA")]
      data_table(data)
    })
    
    #' Generate radiobuttons for the various categories in the pData
    output$colorPoints <- renderUI({
      validate(
        need(input$gene != "", FALSE)
      )
      data <- rmNA(filter_plot_Data())
      # colnames <- names(data)[!names(data) %in% c("mRNA","Sample","status","survival")]
      colnames <- names(data)[!sapply(data, is.numeric)]
      colnames <- colnames[!colnames %in% "Sample"]
      # Create the selectInput for the different pData categories
      div(class="row",
          div(class="col-xs-6",
              selectInput("colorP", "Color by:", choices = c("None",colnames), selected = plot_Type())
          ),
          div(class="col-xs-6",
              selectInput("shapeP", "Shape by:", choices = c("None",colnames), selected = "None")
          )
      )
    })
    
    #' Create the selected plot
    output$plot <- renderPlot({
      # To avoid an error when switching datasets in case the colStrip is not available.
        if(input$typePoint){
        data <- rmNA(filter_plot_Data())
        colnames <- names(data)[!names(data) %in% c("mRNA","Sample","status","survival")]
        validate(
          need(input$colorP %in% c("None", colnames), FALSE) %then%
            need(input$shapeP %in% c("None", colnames), FALSE)
          )
        }
      print(box_Plot())
    }, width = function()ifelse(input$tukeyPlot, input$plot_width * 1.5, input$plot_width), height = function()input$plot_height)
    
    #' Data for the statistic
    observe({filter_plot_Data()})
    dataStat <- reactive({
      mRNA <- filter_plot_Data()[ ,"mRNA"]
      group <- filter_plot_Data()[ ,plot_Type()]
      data <- data.frame(mRNA, group)
      data
    })
    
    #' Summary statistic
    output$summary <- renderTable({    
      data <- dataStat()
      stat <- data.frame(data %>%
                           group_by(group) %>%
                           summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data)[1], 2), "%)" ), # prop.table
                                     median = median (mRNA, na.rm=T), mad = mad(mRNA, na.rm=T), mean = mean(mRNA, na.rm=T), 
                                     sd = sd(mRNA, na.rm=T)))
      row.names(stat) <- stat$group
      tot <- data %>%
        summarise(Sample_count = n(),median = median (mRNA, na.rm=T), 
                  mad = mad(mRNA, na.rm=T), mean = mean(mRNA, na.rm=T), sd = sd(mRNA, na.rm=T))
      stat <- stat[,-1]
      stat <- rbind(stat,TOTAL = tot)
      stat 
    }, align='rrrrrr')
    
    #' Tukey post-hoc test, to combine it with the boxplot and to render in a table
    tukey <- reactive({
      data <-  dataStat()
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff, decreasing = TRUE), ]
      tukey
    })
    
    #' Render tukey
    output$tukeyTest <- renderTable({    
      tukey()
    })
    
    #' Pairwise t test
    output$pairwiseTtest <- renderTable({
      data <-  dataStat()
      pttest <- pairwise.t.test(data$mRNA, data$group, na.rm= TRUE, p.adj = "bonferroni", paired = FALSE)
      pttest$p.value
    })
    
    #' Get the selected download file type.
    download_Plot_FileType <- reactive({
      input$downloadPlotFileType  
    })    
    observe({
      plotFileType    <- input$downloadPlotFileType
      plotFileTypePDF <- plotFileType == "pdf"
      plotUnit    <- ifelse(plotFileTypePDF, "inches", "pixels")
      plotUnitDef <- ifelse(plotFileTypePDF, 7, 600)
      plotUnitMin <- ifelse(plotFileTypePDF, 1, 100)
      plotUnitMax <- ifelse(plotFileTypePDF, 12, 2000)
      plotUnitStep <- ifelse(plotFileTypePDF, 0.1, 50)
      
      updateNumericInput(
        session,
        inputId = "downloadPlotHeight",
        label = sprintf("Height (%s)", plotUnit),
        value = plotUnitDef, min = plotUnitMin, max = plotUnitMax, step = plotUnitStep)
      
      updateNumericInput(
        session,
        inputId = "downloadPlotWidth",
        label = sprintf("Width (%s)", plotUnit),
        value = plotUnitDef, min = plotUnitMin, max = plotUnitMax, step = plotUnitStep)
    })
    
    #' Get the download dimensions.
    download_Plot_Height <- reactive({
      input$downloadPlotHeight
    })
    download_Plot_Width <- reactive({
      input$downloadPlotWidth
    })
    
    #' Download the Plot
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_", input$plotTypeSel,  
               ".", download_Plot_FileType())
      },     
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(download_Plot_FileType())
        plotFunction(file, width = download_Plot_Width(), height = download_Plot_Height())
        print(box_Plot())  
        dev.off(which=dev.cur())
      }
    )
    
    #' Reactive expressions for the conditonal panels to work in the right way
    gcimp_Surv <- reactive ({
      if(input$histologySurv == "GBM") {
        gcimpSurv <- input$gcimpSurv
      } else {
        gcimpSurv <- FALSE
      }
      gcimpSurv
    })
    
    primary_Surv <- reactive ({
      if(input$histologySurv == "GBM") {
        primarySurv <- input$primarySurv
      } else {
        primarySurv <- FALSE
      }
      primarySurv
    })
    
    all_Sub_Surv <- reactive ({
      if(input$histologySurv == "GBM" & input$subtypeSurv == "All") {
        allSubSurv <- input$allSubSurv
      } else {
        allSubSurv <- FALSE
      }
      allSubSurv
    })
    
    #' Extract the survival data.
    surv_Data <- reactive({    
      df <- data()
      df <- filter(df, !is.na(df$status))
      if (input$histologySurv != "All"){
        df <- filter(df, Histology == input$histologySurv)
      }
      if (input$histologySurv == "GBM" & input$subtypeSurv != "All") {
        df <- filter (df, Subtype == input$subtypeSurv)
      }
      # exclude G-CIMP is selected
      if (gcimp_Surv()){
        df <- filter (df, CIMP_status != "G-CIMP")
      }
      # select primary sample
      if (primary_Surv() & any(!is.na(df$Recurrence))) {
        df <- filter (df, Recurrence == "Primary")
      }
      df
    })  
    
    #' Create a slider for the manual cutoff of the Kaplan Meier plot
    mRNA_surv <- reactive({
      surv_Need()
      validate(need(input$histologySurv %in% c("All", histo()),""))      
      mRNA <- surv_Data()[ ,"mRNA"]
      mRNA.values <- round(mRNA[!is.na(mRNA)],2)
      # Generate a vector of continuos values, excluding the first an last value
      mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)]) 
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxmRNA <- renderPlot({    
      validate(need(input$mInput, ""))      
      mRNA <- round(mRNA_surv(),2)
      q <- quantile(mRNA)
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|", col=rgb(0, 0, 0, 0.25))
      # Add a red line to show which  is the current cutoff.
      points(x = input$mInput, y = 0, pch = "|", col="red", cex = 2.5)
      points(x = q[2:4], y = rep(0,3), pch = "|", col="blue", cex = 2)
    }, bg = "transparent")
    
    #' Generate the slider for the manual cutoff
    output$numericCutoff <- renderUI({
      sliderInput(inputId = "mInput",label = NULL, min = min(mRNA_surv()), max = max(mRNA_surv()), 
                  value = median(mRNA_surv()), step = 0.05, round = -2)
    })
    
    #' Requirements for all the survival plots
    surv_Need <- reactive({
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
    
    #' Create a Kaplan Meier plot with cutoff based on quantiles or manual selection
    output$survPlot <- renderPlot({     
      surv_Need ()
      validate(need(input$histologySurv %in% c("All", histo()), FALSE))   
      # Use 'try' to suppress a message throwed the first time manual cutoff is selected
      if(all_Sub_Surv()) {
        try({
          par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2.2,.95,0))
          for (i in c("Classical","Mesenchymal","Neural","Proneural")) {
            survivalPlot (surv_Data(), input$gene, group = "GBM", subtype = i,
                          cutoff = input$cutoff, numeric = input$mInput)
          }
        }, silent = TRUE)} else {
          try(survivalPlot (surv_Data(), input$gene, group = input$histologySurv, subtype = input$subtypeSurv,
                            cutoff = input$cutoff, numeric = input$mInput), silent = TRUE)
        }
    }, height = function(){if(!all_Sub_Surv()) {400} else {650}}, width = function(){if(!all_Sub_Surv()) {500} else {850}})
    
    #' Create a table with the data used in Kaplan Meier plot
    output$survDataTable <- renderDataTable({
      data <- surv_Data()[,c("Sample", "Histology", "Recurrence", "Subtype", "CIMP_status", "mRNA",  "survival", "status")]
      names(data)[7:8] <- c("survival_month", "survival_status")
      data <- rmNA(data)
      strat <- get_cutoff(data$mRNA,input$cutoff,input$mInput)
      if (input$cutoff == "quartiles"){
        strat <- factor(strat,labels = c("1st quartile","2nd quartile","3rd quartile","4th quartile"))
      }
      data <- data.frame(data, cutoff_group = strat)
      data_table(data)
    })
    
    #' Download the survPlot
    output$downloadsurvPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_survPlot.pdf")
      },      
      content = function(file) {
        plotFunction <- match.fun(download_Plot_FileType())
        plotFunction(file, width = download_Plot_Width(), height = download_Plot_Height())
        if(all_Sub_Surv()) {
          par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2.2,.95,0))
          for (i in c("Classical","Mesenchymal","Neural","Proneural")) {
            survivalPlot (surv_Data(), input$gene, group = "GBM", subtype = i,
                          cutoff = input$cutoff, numeric = input$mInput)
          }
        } else {
          survivalPlot (surv_Data(), input$gene, group = input$histologySurv, subtype = input$subtypeSurv,
                        cutoff = input$cutoff, numeric = input$mInput)
        }
        dev.off()
      }
    ) 
    
    #' Subset to GBM samples for the interactive HR plot.
    surv_GBM <- reactive({
      df <- filter (surv_Data(), Histology == "GBM")
    })
    
    #' Extract the GBM expression values for the interactive HR plot.
    gene_Exp <- reactive({
      geneExp <- surv_GBM()[ ,"mRNA"]
      currentClick$stale <<- TRUE
      geneExp
    })
    
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
    
    #' Generate the cutoff value for the interactive HR plot.
    get_Cutoff <- reactive({
      input$hrClick
      gene_Exp()        
      # See if there's been a click since the last gene change.
      if (!is.null(currentClick$click) && !currentClick$stale){
        return(currentClick$click$x)
      }       
      median(gene_Exp())
    })
    
    #' Extract the Hazard ratios for the input gene.
    HR <- reactive ({
      HR <- getHR(surv_GBM())
    })
    
    #' Render a plot to show the the Hazard ratio for the gene's expression values
    output$hazardPlot <- renderPlot({        
      validate(need(!input$dataset %in% c("TCGA LGG","Gorovets","POLA Network"), "Interactive HR plot currently available only for GBM samples") %then%
                 need(histo_Surv_Selected() == "GBM","Please select GBM samples in the 'Histology' dropdown menu") %then%
                 need(!input$dataset %in% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      surv_Need()
      input$tabSurv
      # Plot the hazardplot 
      hazardPlot(HR(), input$quantile)
      # Add a vertical line to show where the current cutoff is.
      abline(v = get_Cutoff(), col = 4)
    }, bg = "transparent")
    
    #' Data used to generate the HR plot
    output$hazardDataTable <- renderDataTable({
      validate(need(!input$dataset %in% c("TCGA LGG","Gorovets","POLA Network"), "Interactive HR plot currently available only for GBM samples") %then%
                 need(histo_Surv_Selected() == "GBM","Please select GBM samples in the 'Histology' dropdown menu") %then%
                 need(!input$dataset %in% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      surv_Need()
      data <- round(HR(),3)
      names(data) <- c("mRNA", "HR", "HR.lower", "HR.upper", "n.obs.1", "n.obs.2")
      data_table(data)
    })
    
    #' A reactive survival formula
    survival_Fml <- reactive({
      # Create the groups based on which samples are above/below the cutoff
      expressionGrp <- as.integer(gene_Exp() < get_Cutoff())
      # Create the survival object 
      surv <- with(surv_GBM(), Surv(survival, status == 1))
      return(surv ~ expressionGrp)
    })
    
    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      validate(
        need(histo_Surv_Selected() == "GBM", FALSE)
      )
      surv_Need()
      cutoff <- get_Cutoff()
      surv <- survival_Fml()
      kmPlot(cutoff, surv)
    })
    
    #' Download the kmPlot
    output$downloadkmPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_kmPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(download_Plot_FileType())
        plotFunction(file, width = download_Plot_Width(), height = download_Plot_Height())
        kmPlot(get_Cutoff(), survival_Fml())
        dev.off()
      }
    )
    
    #' Generate reactive Inputs for the corrPlot to be used also to download the complete plot
    color_by <- reactive({
      switch(input$colorBy, 
             none = "none",
             Histology = "Histology",
             Subtype = "Subtype")
    })
    
    separate_by <- reactive({
      switch(input$separateBy, 
             none = "none",
             Histology = "Histology",
             Subtype = "Subtype")
    })
    
    
    corr_data <- reactive({
      df <- exprs()
      if (input$histologyCorr != "All") {
        df <- filter(df, Histology == input$histologyCorr)
      } 
      if (input$histologyCorr == "GBM" & input$subtype != "All") {
        df <- filter (df, Subtype == input$subtype)
      }
      df
    })
    
    #' Data for the correlation plot
    corr_Two_Genes <- reactive({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        need(input$gene2 != "", "Please enter Gene 2")%then%
          need(input$gene2 %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the choosen histology is not available.
        need(input$histologyCorr %in% c("All",histo()), FALSE)
      )
      data <- corr_data()[,c("Sample", "Histology", "Subtype", input$gene, input$gene2)]
    })
    
    #' Generate the correlation plot
    output$corrPlot <- renderPlot({    
      if (input$dataset == "TCGA LGG") {
        validate(need(input$colorBy != "Subtype" & input$separateBy != "Subtype", "Subtype available for GBM samples only")) 
      }
      p <- myCorggPlot(corr_Two_Genes(), input$gene, input$gene2, colorBy = color_by(), separateBy = separate_by())
      print(p)
    }, width = function()input$cor_plot_width, height = function()input$cor_plot_height)
    
    #' Generate a summary of the correlation test
    output$corrTest <- renderTable({
      df <- corr_Two_Genes()
      if (separate_by() != "none") {
        cor <- substitute(df %>%
                            group_by(g) %>%
                            do(tidy(cor.test(.$x,.$y, use = "complete.obs", method = tolower(input$statCorr)))), 
                          list(g = as.name(input$separateBy), x = as.name(input$gene), y = as.name(input$gene2)))
        cor <- data.frame(eval(cor))
      } else {
        cor <- tidy(cor.test(df[ ,input$gene], df[ ,input$gene2], use = "complete.obs", method = tolower(input$statCorr)))
      }
      cor
    })
    
    #' Table with the data used for the correlation plot
    output$corrDataTable <- renderDataTable({
      data <- corr_Two_Genes()
      data_table(data)
    })
    
    #' Download the corrPlot
    output$downloadcorrPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_corrPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(download_Plot_FileType())
        plotFunction(file, width = download_Plot_Width(), height = download_Plot_Height())
        myCorggPlot (corr_Two_Genes(), input$gene, input$gene2, color_by(), separate_by())
        dev.off()
      }
    )
    
    #' Multiple genes correlation
    corr_Many_Genes <- reactive({
      validate(
        # Need two or more genes
        need(length(input$genelist) > 1, "Please enter 2 or more genes in the panel on the left")%then%
          need(input$genelist %in% names(exprs()),"Gene not available for this dataset")%then%         
          # Trying to avoid an error when switching datasets in case the choosen histology is not available.
          need(input$histologyCorr %in% c("All",histo()), FALSE)
      )
      data <- corr_data()[ ,c("Sample", "Histology", "Subtype", input$genelist)] 
      data
    })
    
    #' Generate the pairs plot
    output$pairsPlot <- renderPlot({
      ggpairs(corr_Many_Genes()[,input$genelist],lower=list(continuous="smooth", params=list(alpha=0.5)))  
    })
    
    #' Table with the data used for the pairs plot
    output$corrPairsDataTable <- renderDataTable({
      data <- corr_Many_Genes()
      data_table(data)
    })
    
    #' Download the pairs plot
    output$downloadpairsPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$dataset, "_pairsPlot.pdf")
      },
      content = function(file) {
        plotFunction <- match.fun(download_Plot_FileType())
        plotFunction(file, width = download_Plot_Width(), height = download_Plot_Height())
        print(ggpairs(corr_Many_Genes()[,input$genelist],
                      lower = list(continuous="smooth", params=c(alpha=0.5)))) 
        dev.off()
      }
    )
    
    #' Generate the correlation table 
    corr_all_data <- reactive ({
      corr <- getCorr(corr_data(), input$gene, tolower(input$corrMethod))
      corr  <- merge(genes, corr, by="Gene")
      corr <- arrange(corr, desc(abs(r)))
    })
    
    #' Generate a reactive element of the the correlation data 
    corr_filter_data <- reactive({  
      corr.table <- suppressWarnings(corr_all_data())  # suppressWarnings  is used to prevent the warning messages in the LGG dataset  
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
    output$corrAllTable <- renderDataTable({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(corr_data()),"Gene not available for this dataset")
      )   
      corr_filter_data()
    }, rownames = FALSE, selection = 'single', extensions = "TableTools", 
    options = list(orderClasses = TRUE, lengthMenu = c(20, 50, 100), pagingType = "full",
                   dom = 'T<"clear">lfrtip', tableTools = list(sSwfPath = copySWF(dest = "www")))
    )
    
    #' Generate a reactive value for the input$rows that set to NULL when the dataset change
    v <- reactiveValues(rows = NULL)
    observeEvent(input$corrAllTable_rows_selected, {
      v$rows <- input$corrAllTable_rows_selected
    })
    observeEvent(datasetInput(), {
      v$rows <- NULL
    })
    observeEvent(input$histologyCorr, {
      v$rows <- NULL
    })
    observeEvent(input$gene, {
      v$rows <- NULL
    })
    
    #' Generate the correlation plot
    output$corrAllPlot <- renderPlot({
      v$rows
      validate(
        need(input$gene != "", FALSE)%then%
          need(v$rows!= "","Click on a row to see the corresponding correlation plot.")
      )
      aes_scatter <- aes_string(x = input$gene, y = v$rows)
      ggplot(corr_data(),mapping = aes_scatter) + theme(legend.position=c(1,1),legend.justification=c(1,1)) +
        geom_point(alpha=.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_linedraw() 
    })
    
    #' RPPA data analysis
    rppa_RNA <- reactive({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG"), FALSE)%then%
          need(input$gene != "", FALSE)%then%
          need(input$gene %in% names(exprs()), FALSE)
      )
      samples <- intersect(row.names(rppas()),exprs()[,"Sample"])
      mRNA <- round(exprs()[samples,input$gene],2)
    })
    
    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxRppaRNA <- renderPlot({    
      validate(
        need(input$rppaCut, FALSE)
      )      
      mRNA <- rppa_RNA()
      q <- quantile(mRNA)
      xrange <-range(mRNA)
      par(mar = c(0,0,0,0)) 
      plot(0, 0, type = "n", xlim = c(xrange[1] + 0.25, xrange[2]) , ylim = c(-0.1,  + 0.1), ylab ="", xlab = "", axes = FALSE)
      points(x = mRNA, y = rep(0, length(mRNA)), pch="|", col=rgb(0, 0, 0, 0.5))
      # Add a red line to show which  is the current cutoff.
      points(x = input$rppaCut, y = 0, pch = "|", col="red", cex = 2.5)
      points(x = q[2:4], y = rep(0,3), pch = "|", col="blue", cex = 2)
      #       abline(v= q[2:4], col="blue")
    }, bg = "transparent")
    
    output$rppaCutoff <- renderUI({
      sliderInput(inputId = "rppaCut",label = "mRNA cutoff", min = min(rppa_RNA()), max = max(rppa_RNA()), 
                  value = median(rppa_RNA()), step = 0.05, round = -2)
    })
    
    output$rppaTable <- renderDataTable({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG"), "RPPA data available only for TCGA GBM and TCGA LGG datasets")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")%then%
          need(input$rppaCut != "", FALSE)%then%
          need(input$rppaCut > sort(rppa_RNA(), FALSE)[2] & input$rppaCut < sort(rppa_RNA(), TRUE)[2], FALSE)
      )
      rppa <- rppas()
      mRNA <- rppa_RNA()
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
      d$adj.p.value <- p.adjust(d$p, method = "bonferroni")
      d <- d[order(d$p.value),]
      d
    }, rownames = FALSE, selection = 'single', options = list(lengthMenu = c(20, 50, 100), pagingType = "full"))
    
    #' Generate a reactive value for the input$rows that set to NULL when the dataset change
    rp <- reactiveValues(rppa.rows = NULL)
    observeEvent(input$rppaTable_rows_selected, {
      rp$rppa.rows <- input$rppaTable_rows_selected
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
    observeEvent(input$rppaCut, {
      rp$rppa.rows <- NULL
    })
    
    #' Generate the RPPA box plot
    output$rppaPlot <- renderPlot({
      input$rppaTable_rows_selected
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG"), FALSE)%then%
          need(input$gene != "", FALSE)%then%
          need(rp$rppa.rows!= "","Click on a row to see the corresponding plots.")%then%
          need(rp$rppa.rows %in% names(rppas()), FALSE)
      )
      mRNA <- rppa_RNA()
      strat <- ifelse(mRNA >= input$rppaCut, c("high"),c("low"))
      strat <- factor(strat,levels = c("low", "high"))
      data <- data.frame(mRNA, strat, rppa = rppas()[,rp$rppa.rows])
      r <- round(cor.test(data$mRNA,data$rppa, use = "complete.obs")$estimate,3)
      p.value <- round(cor.test(data$mRNA,data$rppa, use = "complete.obs")$p.value,4)
      p1 <- ggplot(data, aes(x=strat, y = rppa)) + geom_boxplot(outlier.size = 0) + 
        geom_jitter(aes(colour = strat), position = position_jitter(width = .2), size = 2, alpha = 0.5) + 
        xlab(paste(input$gene, "mRNA")) + ylab(paste(rp$rppa.rows,"RPPA score")) + 
        guides(colour=FALSE) + theme_linedraw()
      p2 <- ggplot(data, aes(x=mRNA, y = rppa)) + geom_point(aes(colour = strat), alpha=.5) +
        geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_linedraw() +
        xlab(sprintf("%s mRNA (log2) \n\n r = %s; p value = %s", input$gene, r, p.value)) + 
        ylab(paste(rp$rppa.rows,"RPPA score")) + theme(legend.position = "none")
      grid.arrange(p1, p2, ncol=1)
    },height = 600)
    
    
    #' Link to cBioportal for mutation analysis
    output$mut_link=renderUI({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG"), "Mutations data available only for TCGA GBM and TCGA LGG datasets")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")
      )
      gene <- input$gene
      if(input$dataset == "TCGA GBM") {
        link <- paste0("http://www.cbioportal.org/index.do?cancer_study_list=gbm_tcga&cancer_study_id=gbm_tcga&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=gbm_tcga_mutations&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=gbm_tcga_sequenced&case_ids=&gene_set_choice=user-defined-list&gene_list=",gene,
                       "%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit")
      } else if (input$dataset == "TCGA LGG") {
        link <- paste0("http://www.cbioportal.org/index.do?cancer_study_list=lgg_tcga&cancer_study_id=lgg_tcga&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=lgg_tcga_mutations&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=lgg_tcga_sequenced&case_ids=&gene_set_choice=user-defined-list&gene_list=",gene,
                       "%0D%0A%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit")
      }
      
      tags$iframe(src= link, style="width: 1000px; height: 600px")
      
      #       ccle <- paste0("http://www.cbioportal.org/index.do?cancer_study_list=cellline_ccle_broad&cancer_study_id=cellline_ccle_broad&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=cellline_ccle_broad_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=cellline_ccle_broad_CNA&genetic_profile_ids_PROFILE_MRNA_EXPRESSION=cellline_ccle_broad_mrna_median_Zscores&Z_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=cellline_ccle_broad_central_nervous_system&case_ids=&gene_set_choice=user-defined-list&gene_list=",gene,
      #                      "&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit")
      #       
      #       cosmic <- paste0("http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=", gene)
      #       
      #       div(
      #         a('TCGA GBM', href = gbm_mut),
      #         br(),
      #         a('TCGA LGG', href = lgg_mut),
      #         br(),
      #         a('Cancer Cell Line Encyclopedia', href = ccle),
      #         br(),
      #         a('COSMIC', href = cosmic)
      #       )
    })
    
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
            t <- tableGrob(stat)
                           # , gp = gpar(fontsize=14),row.just = "right", core.just = "right")
#             t <- tableGrob(stat, theme = ttheme_default(core = list(fg_params = list(hjust=1, x=0.9)),
#                                                         rowhead = list(fg_params=list(hjust=1, x=0.95))))
            grid.arrange(p, t, ncol = 2, just = c("center", "top")) # `just` it's not working
          })
        })
      }
    })
    
    #' Reactive function for an HTML table view of the data
    summary_data <- reactive({
      mRNA <- exprs()[ , input$gene, drop = FALSE]
      names(mRNA) <- paste0(input$gene,"_mRNA")
      data <- cbind(pDatas(),mRNA)
      if (input$dataset %in% c("TCGA GBM", "TCGA LGG", "TCGA GBMLGG")) {
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
        data <- summary_data()
      }
      datatable(data, rownames = FALSE, extensions = c("FixedColumns", "TableTools"),
                options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE,
                               lengthMenu = list(c(10, 25, 50, -1), c('10','25','50','All')), dom = 'T<"clear">lfrtip',
                               tableTools = list(sSwfPath = copySWF(dest = "www"))))
    })
    
    #' Generate a graphic summary of the dataset, using ggvis
    output$piePlots <- renderUI({
      data <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status")]
      data <- rmNA(data) # Removing unavailable (all NA) groups
      plot_output_list <- lapply(names(data), function(i) {
        plotname <- paste("plot", i, sep="")
        htmlOutput(plotname)
      })
      # Convert the list to a tagList - this is necessary for the list of items to display properly.
      do.call(tagList, plot_output_list)
    })
    
    observe ({                                                               
      data <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status")]
      data <- rmNA(data)
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
      df <-  rmNA(df) # Removing unavailable (all NA) groups
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
        Training <- train$Subtype # Expression subtype
        subtypes <- c("Classical","Mesenchymal","Neural","Proneural")
      } else if (input$tumorType == "lgg") {
        train <- lgg.tcga[["expr"]]
        train <- train[train$Sample %in% lgg.core.samples,]
        pData <- lgg.tcga[["pData"]]
        pData.core <- pData[row.names(train),]
        Training <- pData.core$Subtype
        subtypes <- c("IDHmut-codel","IDHmut-non-codel","IDHwt")
        #         Training <- pData.core$RNASeqCluster # Expression subtype (not the molecular subtype)
        #         subtypes <- c("R1","R2","R3","R4")
      }
      row.names(train) <- train[,"Sample"]
      train.exp <- data.frame(t(train[,-c(1:8)])) #select only the expression data
      row.names(upData) <- upData[,"Sample"]
      upData <- upData [,-1] #remove the column with Sample ID
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
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goSvm != 0,'Please press "Submit SVM"')
      )
      data_table(svm.call())
    })
    
    #' Reactive function to generate k-nearest neighbour subtype call to pass to data table and download handler
    knn.call <- eventReactive (input$goKnn | input$goSub3, {
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      if (input$tumorType == "gbm") {
        train <- gbm.tcga[["expr"]]
        train <- train[train$Sample %in% gbm.core.samples,]
        Training <- train$Subtype # Expression subtype
        k <- 15
        subtypes <- c("Classical","Mesenchymal","Neural","Proneural")
      } else if (input$tumorType == "lgg") {
        train <- lgg.tcga[["expr"]]
        train <- train[train$Sample %in% lgg.core.samples,]
        pData <- lgg.tcga[["pData"]]
        pData.core <- pData[row.names(train),]
        Training <- pData.core$Subtype
        subtypes <- c("IDHmut-codel","IDHmut-non-codel","IDHwt")
#         Training <- pData.core$RNASeqCluster # Expression subtype (not the molecular subtype)
#         subtypes <- c("R1","R2","R3","R4")
        k <- 3
      }
      row.names(train) <- train[,"Sample"]
      train.exp <- as.matrix(train[,-c(1:8)])
      row.names(upData) <- upData[,"Sample"]
      learn.exp <-as.matrix(upData [,-1])
      # Common genes of the two datasets
      genes <- intersect(colnames(train.exp), colnames(learn.exp))
      set.seed(1234)
      pred <- knn3Train(train = train.exp[,genes], test = learn.exp[,genes], cl =  Training, k = k, prob=TRUE)
      kn <- data.frame(Sample = rownames(upData), knn.subtype.call = pred, prob = round(attr(pred,"prob"),2))
      names(kn)[3:5] <- subtypes
      kn
    })
    
    #' Rerndering the knn subtype call as a data table
    output$knn <- renderDataTable({
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goKnn != 0,'Please press "Submit K-NN"')
      )
      data_table(knn.call())
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
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq=FALSE, parallel.sz = 0,
                           min.sz=0, max.sz=10000, verbose=TRUE)
      subtype_scores <- round(t(gsva_results),3)
      subtype_final <- data.frame(Sample = rownames(upData), gsea.subtype.call = names(gene_list)[apply(subtype_scores,1,which.max)], 
                                  subtype_scores)
      subtype_final <- subtype_final[,c("Sample", "gsea.subtype.call", "Classical","Mesenchymal","Neural","Proneural")]
      subtype_final
    })
    
    #' Rerndering the subtype call as a data table
    output$gsva <- renderDataTable({ 
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goGsva != 0,'Please press "Submit ssGSEA"')
      )
      data_table(gsva.call())
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
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goSub3 != 0, 'Please press "Submit 3-Way"')
      )
      data_table(sub3.call())
    })
    
    output$sub3Plot <- renderPlot({
      data <- sub3.call()[,1:4]
      data$Sample <- factor(data$Sample, levels=(data$Sample)[order(data$svm.call)]) # sort the sample on the svm.call
      data <- melt(data, id.var = "Sample")
      names(data) <- c("Sample","Call","Subtype")
      plot <- ggplot(data, aes(Sample, Call)) + geom_tile(aes(fill = Subtype), colour = "white") + ylab("") + xlab("") +
        scale_x_discrete(expand = c(0, 0)) + theme_minimal() + scale_fill_brewer(palette = "Set1") + # scale_fill_manual(values = terrain.colors(4)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(legend.position="bottom")
      #       table <- tableGrob(table(data$Call,data$Subtype), gp = gpar(fontsize=11),
      #                      row.just = "right", core.just = "right")
      #       h <- grobHeight(table)
      #       w <- grobWidth(table)
      #       title <- textGrob("Subtype call summary:", x=unit(0.5,"npc") - 0.3*w, 
      #                         y=unit(0.5,"npc") + 0.6*h, vjust=0, gp=gpar(fontsize=12))
      #       gt <- gTree(children=gList(table, title))
      #       grid.arrange(plot1, gt, ncol = 2, widths = c(6,2))
      plot
    }, height = 200)
    
    #' Reactive function to generate Estimate call to pass to data table 
    est.call <- eventReactive (input$goEst,{
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      row.names(upData) <- upData[,"Sample"]
      ds <- data.frame(t(upData[,-1]))
      est <- myEstimateScore(ds, platform = input$platformEst)
      est <- cbind(Sample = upData$Sample, round(est,2))
    })
    
    #' Rerndering the subtype call as a data table
    output$estScore <- renderDataTable({
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goEst != 0,'Please press "Submit ESTIMATE"')
      )
      datatable(est.call(), rownames = FALSE, selection = 'single', extensions = "TableTools", 
                options = list(orderClasses = TRUE, lengthMenu = list(c(20, 50, 100, -1), c('20','50','100','All')), pageLength = 20, pagingType = "full", autoWidth = TRUE,
                               dom = 'T<"clear">lfrtip', tableTools = list(sSwfPath = copySWF(dest = "www"))))
    })
    
    #' Generate the Purity plot
    output$purityPlot <- renderPlot({
      validate(
        need(est.call(), FALSE)%then%
          need(input$platformEst == "affymetrix", "Sorry, the purity plot for Agilent and RNAseq have not yet been impemented") %then%
          need(input$estScore_rows_selected != "","Click on a row to see the corresponding purity plot.")       
      )
      plotPurity(est.call(), input$estScore_rows_selected, platform = input$platformEst)
    })
       
    #' Reactive function to generate Deconvolute scores to pass to data table 
    deconv.call <- eventReactive (input$goDec,{
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      row.names(upData) <- upData[,"Sample"]
      exprs <- data.frame(t(upData[,-1]),check.names=FALSE)
      if (input$geneListDec == "engler") {
        gene_list <- engler_gene_set_list
      } else if (input$geneListDec == "galon") {
        gene_list <- galon_gene_set_list
      } else if (input$geneListDec == "galon_engler") {
        gene_list <- galon_engler
      }
      platformDec <- ifelse(input$platformDec == "rnaseq",TRUE,FALSE)
      set.seed(1234)
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq = platformDec, parallel.sz = 0,
                           min.sz=25, max.sz=10000, verbose=TRUE)
      deconv_scores <- data.frame(Sample = rownames(upData),t(gsva_results))
      deconv <- list(results = gsva_results, scores = deconv_scores)
    })
    
    output$pDataDec <- renderUI({
      validate(
        need(!is.null(input$pDataFile),"Please upload the pData to be included")
      )
      inFile <- input$pDataFile
      pData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      pData <- rmNA(pData)
      pData_group <- names(pData)
      pData_group <- pData_group[pData_group!="Sample"]
      # Create the selectInput for the different pData categories
      selectInput(inputId = "pDataDec", label = "Select pData group", choices = pData_group, selectize = TRUE)
    })
    
    #' Required for the conditional panel 'deconvPData' to work correctly
    observe({
      if(input$tabTools != "DeconvoluteMe")
        updateCheckboxInput(session, inputId = "deconvPData", value = FALSE) 
    })
    
    #' Rerndering the Deconvolute scores as a heatmap
    output$deconvHeatmap <- renderPlot({ 
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goDec != 0,'Please press "Submit Deconvolute"')
      )
      data <- deconv.call()[["results"]]
      heatmap3(data,margins = c(0,11),labCol=NA,cexRow = 1.25)
      if(input$deconvPData){
        validate(
          need(!is.null(input$pDataFile),"Please upload the pData to be included")
        )
        inFile <- input$pDataFile
        pData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
        pData$Sample <- as.character(pData$Sample)
        validate(
          need(all.equal(pData$Sample,names(data)),"ERROR: Samples don't match")
        )
        pData <- rmNA(pData)
        n <- nlevels(pData[,input$pDataDec])
        colors <- data.frame(levels(pData[,input$pDataDec]), color = brewer_pal(palette = "Dark2")(n))
        names(colors)[1] <- paste0(input$pDataDec)
        dataColor <- merge(pData, colors, all=T)
        dataColor <- dataColor[order(dataColor$Sample),]
        col <- as.character(dataColor$color)
        heatmap3(data,margins = c(0,11),labCol=NA,cexRow = 1.25,ColSideColors = col,ColSideLabs = paste0(input$pDataDec))
        legend(0.925, 1, legend = colors[,input$pDataDec], fill = as.character(colors[,"color"]),
               border = FALSE, bty = "n", y.intersp = 1, cex = 0.8, title = paste0(input$pDataDec))
      }
    })
    
    #' Rerndering the Deconvolute scores as boxplot
    output$deconvBoxPlot <- renderPlot({ 
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goDec != 0,'Please press "Submit Deconvolute"')
      )
      data <- deconv.call()[["results"]]
      melted <- melt(data)
      p <- ggplot(melted,aes(x=Var1,y=value,fill=Var1)) + geom_boxplot() + xlab("") + ylab("score") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.title = element_blank())
      if(input$deconvPData){
        validate(
          need(!is.null(input$pDataFile),"Please upload the pData to be included")
        )
        inFile <- input$pDataFile
        pData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
        pData$Sample <- as.character(pData$Sample)
        validate(
          need(all.equal(pData$Sample,names(data)),"ERROR: Samples don't match")
        )
        melted <- melt(merge(t(data),pData[,c("Sample",input$pDataDec)],
                             by.x="row.names",by.y="Sample"))
        p <- ggplot(melted,aes_string(x=input$pDataDec,y="value",fill=input$pDataDec)) + geom_boxplot() + xlab("") + ylab("score") +
          facet_wrap(~variable)
      }
      p
    },height = 800)
    
    #' Rerndering the Deconvolute scores as a data table
    output$deconvScore <- renderDataTable({ 
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goDec != 0,'Please press "Submit Deconvolute"')
      )
      datatable(deconv.call()[["scores"]], rownames = FALSE, extensions = c("FixedColumns", "TableTools"),
                options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE,
                               lengthMenu = list(c(20, 50, 100, -1), c('20','50','100','All')), 
                               pageLength = 20, dom = 'T<"clear">lfrtip',
                               tableTools = list(sSwfPath = copySWF(dest = "www")))
                )
    })
    
  })
