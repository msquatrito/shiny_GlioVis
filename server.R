#' server.R for Gliovis

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
`%notin%` <- Negate('%in%')
shinyServer(
  function(input, output, session) {

    options(shiny.maxRequestSize=100*1024^2)

    dataset_type <- reactive({
      switch(input$datasetType,
             "Adult" = adult_datasets,
             "Pediatric" = pediatric_datasets)
    })

    observe({
      updateSelectizeInput(session, inputId = "dataset", choices = dataset_type(), selected = dataset_type()[1], server = TRUE)
    })

    #' Return the requested dataset
    dataset_input <- reactive({
      switch(input$dataset,
             "TCGA GBM" = gbm.tcga,
             "TCGA LGG" = lgg.tcga,
             "TCGA GBMLGG" = lgg_gbm.tcga,
             "Rembrandt" = rembrandt,
             "Gravendeel" = gravendeel,
             "Kamoun" = kamoun,
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
             "Walsh" = walsh,
             "Grzmil" = grzmil,
             "Donson" = donson,
             "Li" = li,
             "Vital" = vital,
             "Joo" = joo,
             "Oh" = oh,
             "Ivy GAP" = ivy,
             "POLA Network" = pola,
             "Kwom" = kwom,
             "Gleize" = gleize,
             "Sturm_2016" = sturm_2016,
             "Paugh" = paugh,
             "Mascelli" = mascelli,
             "Schwartzentruber" = schwartzentruber,
             "Lambert" = lambert,
             "Griesinger" = griesinger,
             "Zakrzewski" = zakrzewski,
             "Sturm_2012" = sturm_2012,
             "Bender" = bender,
             "Bergthold" = bergthold,
             "Henriquez" = henriquez,
             "Buczkowicz" = buczkowicz,
             "Witt" = witt,
             "deBont" = deBont,
             "Gump" = gump,
             "Hoffman" = hoffman,
             "Johnson" = johnson,
             "Northcott_2012" = northcott_2012,
             "Northcott_2011" = northcott_2011,
             "Pomeroy" = pomeroy,
             "Kool" = kool,
             "Johann" = johann)
             
    })

    #' Expression data
    exprs <- reactive({
      if(input$dataset == "TCGA GBM") {
        switch(input$microarray_platform,
        "HG-U133A" = dataset_input()[["expr"]],
        "Agilent-4502A" =dataset_input()[["agilent"]],
        "RNA-Seq" = dataset_input()[["rseq"]]
        )
      } else {
        dataset_input()[["expr"]]
      }
    })

    #' Phenotype data
    pDatas <- reactive({
      dataset_input()[["pData"]]
    })

    #' CNA data
    cnas <- reactive({
      dataset_input()[["cna"]]
    })

    #' RPPA data
    rppas <- reactive({
      dataset_input()[["rppa"]]
    })

    #' Text matching with the gene names list
    updateSelectizeInput(session, inputId = "gene", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "gene2", choices = gene_names, server = TRUE)
    updateSelectizeInput(session, inputId = "corrGene", choices = gene_names, server = TRUE)
    observe({
      updateSelectizeInput(session, inputId = "mutGene", choices = gene_names, selected = input$gene, server = TRUE)
    })


    #' Required for the conditional panel 'corrMany' to work correctly
    observe({
      if(input$tab1 != 3)
        updateTabsetPanel(session, inputId = "tabCorr", selected = "corrTwo")
    })

    observe({
      if(input$tabTools != "DeconvoluteMe")
        updateCheckboxInput(session, inputId = "deconvPData", value = FALSE)
    })

    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plot_selected <- reactive ({
      if (input$plotTypeSel %in% plotList[[input$dataset]]){
        input$plotTypeSel
      } else {
        NULL
      }
    })

    #' Return the names of the available user-defined plots
    plot_user_selection <- reactive ({
      # Exclude the pre-defeined plots and numeric variabe
      dropit <- c("Sample","Histology","Grade","Recurrence","Subtype", "CIMP_status", "survival",
                  "status", "Age", "ID","Patient_ID","Sample_ID", "Matching.sample", "Therapy_Class","title","tumor_name")
      data <- pDatas()[,names(pDatas()) %notin% dropit, drop = FALSE]
      n <- names(data)
      n
    })

    #' When switching datasets, if the selected plot is not available it will choose the first plot of the list
    plot_user_selected <- reactive({
      if (input$plotTypeUserSel %in% plot_user_selection()){
        input$plotTypeUserSel
      } else {
        NULL
      }
    })

    plot_type <- reactive({
      if (input$plotType == "Pre-defined"){
        plot_selected()
      } else if (input$plotType == "User-defined") {
        plot_user_selected()
      }
    })

    #' Change the plot type available for a specific dataset
    observe({
      updateSelectInput(session, inputId = "plotTypeSel", choices = plotList[[input$dataset]], selected = plot_selected())
      updateSelectInput(session, inputId = "plotTypeUserSel", choices = plot_user_selection(), selected = plot_user_selected())
    }, priority = 10)

    # Caption with gene and dataset
    output$caption <- renderText({
      if (input$gene == "" )
        return()
      title <- paste(input$gene, "in", input$dataset, "dataset")
    })

    #' Return the available histology, to be used in the updateSelectInput
    histo <- reactive({
      levels(exprs()[,"Histology"])
    })

    #' Return the available subtype, to be used in the updateSelectInput
    subtype <- reactive({
      req(input$dataset)
      if (input$histology == "All"){
        df <- exprs()
        subtype <- levels(df$Subtype)
      } else{
        df <- subset(exprs(), Histology == input$histology)
        subtype <- intersect(levels(df$Subtype),df$Subtype)
      }
      subtype
    })

    #' When switching datasets if the selected histo is not available it will choose "All"
    histo_selected <- reactive ({
      if (input$histology %in% c("All", histo())){
        input$histology
      } else {
        "All"
      }
    })

    #' When switching datasets if the selected subtype is not available it will choose "All"
    subtype_selected <- reactive ({
      if (input$subtype %in% c("All", subtype())){
        input$subtype
      } else {
        "All"
      }
    })

    observe({
      updateSelectInput(session, inputId = "histology", choices = c("All", histo()), selected = histo_selected())
      updateSelectInput(session, inputId = "subtype", choices = c("All", subtype()), selected = subtype_selected())
    })

    #' Generate a dataframe with the data to plot
    data <- reactive({
      # Trying to avoid an error when switching datasets in case the plotType is not available.
      req(c(plot_selected(),plot_user_selected()) %in% c(plotList[[input$dataset]],plot_user_selection()))
      # req(plot_user_selected() %in% plot_user_selection())
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")
      )
      mRNA <- exprs()[ ,input$gene]
      if (input$scale) {
        mRNA <- scale(mRNA)
      }
      data <- cbind(mRNA, exprs()[,2:6]) # To combine with pData
      data <- cbind(data, pDatas()[,!names(pDatas()) %in% names(data)]) # To combine with more pData for the report
      if (input$dataset %in% c("TCGA GBM", "TCGA LGG", "TCGA GBMLGG")) {
        if(input$gene %in% names(cnas())){
          Copy_number <- cnas()[ ,input$gene]
          Copy_number <- factor(Copy_number, levels = c(-2:2), labels = c("Homdel", "Hetloss", "Diploid", "Gain", "Amp"))
          Copy_number <- droplevels(Copy_number)
        } else {
          Copy_number <- NA # Some genes don't have copy numbers data
        }
        data <- cbind(Copy_number,data)
      }
      data <- data[complete.cases(data[,"mRNA"]),] # Aiglent TCGA GBM data doesn't have all the observations
      data
    })

    #' Data for the box plot
    plot_data <- reactive({
      validate(
        need(!all(is.na(data()[ ,plot_type()])),"Sorry,no gene data available for this group")
      )
      data <- data()
      data <- subset(data, !is.na(data[ ,plot_type()]))
    })

    #' Filtered data for the box plot
    filter_plot_data <- reactive({
      if (input$removeMe) {
        req(input$removeGp %in% levels(plot_data()[ ,plot_type()])) #To silence an error thrown by the stat analysis
        data <- plot_data()
        data <- subset(data, data[ ,plot_type()] %in% input$removeGp)
        data[ ,plot_type()] <- factor(data[ ,plot_type()],levels = input$removeGp)
      } else {
        data <- plot_data()
      }
      data
    })

    #' Select input for groups to add/exclude
    output$removeGp <- renderUI({
      groups <- levels(plot_data()[ ,plot_type()])
      selectizeInput(inputId = "removeGp", label = "", choices = groups, selected = groups, multiple = TRUE,
                     options = list(plugins = list('remove_button','drag_drop')))
    })

    #' Reset input$removeMe when plot_data() changes
    observeEvent(plot_data(), {
      updateCheckboxInput(session, inputId = "removeMe", value = FALSE)
    })

    #' Populate Xaxis labels
    observe({
      updateTextInput(session, inputId = "myXlab",value = paste0("\n",plot_type()))
    })

    # Tukey plot active only when tukey stat data are shown
    observeEvent(!input$tukeyHSD, {
      updateCheckboxInput(session, inputId = "tukeyPlot", value = FALSE)
    })

    #' Reactive function to generate the box plots
    box_Plot <- reactive ({
      data <- filter_plot_data()
      xlabel <- paste("\n", input$myXlab)
      ylabel <- paste(input$myYlab,"\n")
      fillBox <- NULL
      col <- input$colorP
      shape <- input$shapeP
      theme <- theme(axis.text.x = element_text(size = input$axis_text_size, angle = input$xaxisLabelAngle, hjust = ifelse(input$xaxisLabelAngle == 0,0.5,1)),
                     axis.text.y = element_text(size = input$axis_text_size),
                     axis.title.x = element_text(size = input$axis_title_size),
                     axis.title.y = element_text(size = input$axis_title_size))

      if (input$colBox) fillBox <- plot_type()
      if(input$colorP == "None") col <-  NULL
      if(input$shapeP == "None") shape <-  NULL
      if (input$bw) theme <- theme_bw () + theme

      p <- ggplot(data, mapping=aes_string(x=plot_type(), y = "mRNA")) +
        geom_boxplot(aes_string(fill = fillBox), outlier.size = 0, outlier.stroke = 0) +
        geom_jitter(position = position_jitter(width = .5), mapping = aes_string(colour = col, shape = shape),
                    size = input$point_size, alpha = input$alpha) +
        ylab(ylabel) + xlab(xlabel) + theme

      if (input$tukeyPlot) {
        t <- tukey() %>%
          mutate(comparison = row.names(.)) %>%
          ggplot(aes(reorder(comparison, diff), diff, ymin = lwr, ymax= upr, colour = Significance)) +
          geom_point() + geom_errorbar(width = 0.25) + ylab("\nDifferences in mean levels") + xlab("") +
          geom_hline(yintercept = 0, colour="darkgray", linetype = "longdash") + coord_flip() + theme
        p <- grid.arrange(p, t, ncol=2, widths = c(3,2))
      }

      return(p)

    })

    #' Table with the data used for the plot
    output$filterDataTable <- renderDataTable({
      data <- filter_plot_data()[,c("Sample", plot_type(), "mRNA")]
      data_table(data)
    })

    observe({
      data <- rmNA(filter_plot_data())
      colnames <- names(data)[!sapply(data, is.numeric)] # remove muneric categories
      colnames <- colnames[colnames %notin% "Sample"]
      # Create the selectInput for the different pData categories
      updateSelectInput(session, "colorP", "Color by:", choices = c("None",colnames), selected = "None")
      updateSelectInput(session, "shapeP", "Shape by:", choices = c("None",colnames), selected = "None")
    })

    box_width <- reactive({
      if(input$tukeyPlot)
        input$plot_width* 2 else
          input$plot_width
    })

    #' Create the selected plot
    output$plot <- renderPlot({
      # To avoid an error when switching datasets in case the colStrip is not available.
      data <- rmNA(filter_plot_data())
      colnames <- names(data)[names(data) %notin% c("mRNA","Sample","status","survival")]
      req(input$colorP %in% c("None", colnames) & input$shapeP %in% c("None", colnames))
      box_Plot()
    }, width = box_width, height = function()input$plot_height)

    #' Data for the statistic
    observe({filter_plot_data()})
    stat_data <- reactive({
      mRNA <- filter_plot_data()[ ,"mRNA"]
      group <- filter_plot_data()[ ,plot_type()]
      data <- data.frame(mRNA, group)
      data
    })

    #' Summary statistic
    output$summary <- renderTable({
      data <- stat_data()
      stat <- data %>%
        group_by(group) %>%
        summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data)[1], 2), "%)" ), # prop.table
                  median = median(mRNA, na.rm=T), mad = mad(mRNA, na.rm=T), mean = mean(mRNA, na.rm=T),
                  sd = sd(mRNA, na.rm=T)) %>%
        data.frame()
      row.names(stat) <- stat$group
      tot <- data %>%
        summarise(Sample_count = n(), median = median(mRNA, na.rm=T),
                  mad = mad(mRNA, na.rm=T), mean = mean(mRNA, na.rm=T), sd = sd(mRNA, na.rm=T))
      stat <- stat[,-1]
      stat <- rbind(stat,TOTAL = tot)
      stat
    }, align='rrrrrr')

    #' Tukey post-hoc test, to combine it with the boxplot and to render in a table
    tukey <- reactive({
      validate(
        need(nlevels(stat_data()$group)>1,message = "There is only one category, group comparison cannot be performed")
      )
      data <-  stat_data()
      tukey <- data.frame(TukeyHSD(aov(mRNA ~ group, data = data))[[1]])
      tukey$Significance <- as.factor(starmaker(tukey$p.adj, p.levels = c(.001, .01, .05, 1), symbols=c("***", "**", "*", "ns")))
      tukey <- tukey[order(tukey$diff, decreasing = TRUE), ]
      tukey
    })

    #' Render tukey
    output$tukeyTest <- renderTable({
      tukey()
    }, digits = c(2,2,2,2,-1,2))

    #' Pairwise t test
    output$pairwiseTtest <- renderTable({
      validate(
        need(nlevels(stat_data()$group)>1,message = "There is only one category, group comparison cannot be performed")
      )
      data <-  stat_data()
      pttest <- pairwise.t.test(data$mRNA, data$group, na.rm= TRUE, p.adj = "bonferroni", paired = FALSE)
      pttest$p.value
    }, digits = -1)

    #' Get the selected download file type.
    download_plot_file_type <- reactive({
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
    download_plot_height <- reactive({
      input$downloadPlotHeight
    })
    download_plot_width <- reactive({
      input$downloadPlotWidth
    })

    #' Download the Plot
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_", input$plotTypeSel,
               ".", download_plot_file_type())
      },
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height =  download_plot_height())
        if (input$tukeyPlot) {
          grid.draw(box_Plot())
        } else {
          print(box_Plot())
        }
        dev.off()
      }
    )

    #' Observers for the conditonal panels to work in the proper way
    observeEvent(input$histology != "GBM", {
      updateCheckboxInput(session, "primarySurv", value = FALSE)
      updateCheckboxInput(session, "gcimpSurv", value = FALSE)
    })
    observeEvent(input$dataset =="TCGA GBMLGG", {
      updateCheckboxInput(session, "gcimpSurv", value = FALSE)
    })
    observeEvent(input$subtype != "All", {
      updateCheckboxInput(session, "allSubSurv", value = FALSE)
    })

    #' Extract the survival data.
    surv_data <- reactive({
      df <- data()
      df <- subset(df, !is.na(df$status))
      if (input$histology != "All"){
        df <- subset(df, Histology == input$histology)
      }
      if (input$subtype != "All") {
        df <- subset(df, Subtype == input$subtype)
      }
      # exclude G-CIMP is selected
      if (input$gcimpSurv){
        df <- subset(df, CIMP_status != "G-CIMP")
      }
      # select primary sample
      if (input$primarySurv & any(!is.na(df$Recurrence))) {
        df <- subset(df, Recurrence == "Primary")
      }
      data <- df[,c("Sample", "Histology", "Recurrence", "Subtype", "CIMP_status", "mRNA",  "survival", "status")]
      names(data)[7:8] <- c("survival_month", "survival_status")
      data <- rmNA(data)
    })

    #' Create a slider for the manual cutoff of the Kaplan Meier plot
    mRNA_surv <- reactive({
      surv_need()
      req(input$histology %in% c("All", histo()))
      mRNA <- surv_data()[ ,"mRNA"]
      mRNA.values <- round(mRNA[!is.na(mRNA)],2)
      # Generate a vector of continuos values, excluding the first an last value
      mRNA.values <- sort(mRNA.values[mRNA.values != min(mRNA.values) & mRNA.values != max(mRNA.values)])
    })

    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxmRNA <- renderPlot({
      req(input$mInput)
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
    surv_need <- reactive({
      validate(
        need(input$dataset %notin% no_surv_dataset, "Sorry, no survival data are available for this dataset")%then%
          need(input$histology != "Non-tumor","Sorry, no survival data are available for this group")%then%
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

    #' Set survival plot height
    surv_plot_height <- reactive({
      if(input$allSubSurv){
        ifelse(length(subtype())>4, 1300, 650)
      } else {
        ifelse(input$riskTable, 500,325)
      }
    })

    observeEvent(input$allSubSurv, {
      updateCheckboxInput(session, inputId = "riskTable", value = FALSE)
    })

    #' Create a Kaplan Meier plot with cutoff based on quantiles or manual selection
    output$survPlot <- renderPlot({
      surv_need ()
      req(input$histology %in% c("All", histo()))
      # Use 'try' to suppress a message throwed the first time manual cutoff is selected
      if(input$allSubSurv) {
        nrow <- ceiling(length(subtype())/2)
        try({
          p <- list()
          for (i in subtype()) {
            p[[i]] <- survivalPlot(surv_data(), input$gene, group = input$histology, subtype = i, font.legend = input$surv_legend_size,
                                   input$riskTable, cutoff = input$cutoff, numeric = input$mInput, censor = input$censor, conf.int = input$confInt)$plot
          }
          do.call(grid.arrange,args = list(grobs = p, nrow = nrow, ncol=2))
        }, silent = TRUE)
      } else {
        try(survivalPlot(surv_data(), input$gene, group = input$histology, subtype = input$subtype, font.legend = input$surv_legend_size,
                         input$riskTable, cutoff = input$cutoff, numeric = input$mInput, censor = input$censor, conf.int = input$confInt), silent = TRUE)
      }
    }, height = surv_plot_height, width = function(){if(!input$allSubSurv) {500} else {900}})


    #' Create a table with the data used in Kaplan Meier plot
    output$survDataTable <- renderDataTable({
      data <- surv_data()
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
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_survPlot.", download_plot_file_type())
      },
      content = function(file) {
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height =  download_plot_height())
        if(input$allSubSurv) {
          nrow <- ceiling(length(subtype())/2)
            p <- list()
            for (i in subtype()) {
              p[[i]] <- survivalPlot(surv_data(), input$gene, group = input$histology, subtype = i, font.legend = input$surv_legend_size,
                                     input$riskTable, cutoff = input$cutoff, numeric = input$mInput, censor = input$censor, conf.int = input$confInt)$plot
            }
            do.call(grid.arrange,args = list(grobs = p, nrow = nrow, ncol=2))
        } else {
          survivalPlot(surv_data(), input$gene, group = input$histology, subtype = input$subtype, font.legend = input$surv_legend_size,
                        input$riskTable, cutoff = input$cutoff, numeric = input$mInput, censor = input$censor, conf.int = input$confInt)
        }
        dev.off()
      }
    )

    #' Subset to GBM samples for the interactive HR plot.
    surv_GBM <- reactive({
      df <- filter(surv_data(), Histology == "GBM")
    })

    #' Extract the GBM expression values for the interactive HR plot.
    gene_exp <- reactive({
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
      gene_exp()
      # See if there's been a click since the last gene change.
      if (!is.null(currentClick$click) && !currentClick$stale){
        return(currentClick$click$x)
      }
      median(gene_exp())
    })

    #' Extract the Hazard ratios for the input gene.
    HR <- reactive ({
      HR <- getHR(surv_GBM())
    })

    #' Render a plot to show the the Hazard ratio for the gene's expression values
    output$hazardPlot <- renderPlot({
      validate(need(input$dataset %notin% c("TCGA LGG","Gorovets","POLA Network"), "Interactive HR plot currently available only for GBM samples") %then%
                 need(histo_selected() == "GBM","Please select GBM samples in the 'Histology' dropdown menu") %then%
                 need(input$dataset %notin% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      surv_need()
      input$tabSurv
      # Plot the hazardplot
      hazardPlot(HR(), input$quantile)
      # Add a vertical line to show where the current cutoff is.
      abline(v = get_Cutoff(), col = 4)
    }, bg = "transparent")

    #' Data used to generate the HR plot
    output$hazardDataTable <- renderDataTable({
      validate(need(input$dataset %notin% c("TCGA LGG","Gorovets","POLA Network"), "Interactive HR plot currently available only for GBM samples") %then%
                 need(histo_selected() == "GBM","Please select GBM samples in the 'Histology' dropdown menu") %then%
                 need(input$dataset %notin% c("Grzmil","Vital"), "Sorry, too few samples to properly render the HR plot"))
      surv_need()
      data <- round(HR(),3)
      names(data) <- c("mRNA", "HR", "HR.lower", "HR.upper", "n.obs.1", "n.obs.2")
      data_table(data)
    })

    #' A reactive survival formula
    survival_Fml <- reactive({
      # Create the groups based on which samples are above/below the cutoff
      expressionGrp <- as.integer(gene_exp() < get_Cutoff())
      # Create the survival object
      surv <- with(surv_GBM(), Surv(survival_month, survival_status == 1))
      return(surv ~ expressionGrp)
    })

    #' Create a Kaplan Meier plot on the HR cutoff
    output$kmPlot <- renderPlot({
      req(histo_selected() == "GBM")
      surv_need()
      kmPlot(cutoff = get_Cutoff(), surv = survival_Fml(), censor = input$censor, conf.int = input$confInt, font.legend = input$surv_legend_size, input$riskTable)
    }, height = surv_plot_height)

    #' Download the kmPlot
    output$downloadkmPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_kmPlot.",download_plot_file_type())
      },
      content = function(file) {
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height =  download_plot_height())
        kmPlot(cutoff = get_Cutoff(), surv = survival_Fml(), censor = input$censor, conf.int = input$confInt, font.legend = input$surv_legend_size, input$riskTable)
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
      if (input$histology != "All") {
        df <- subset(df, Histology == input$histology)
      }
      if (input$subtype != "All") {
        df <- subset(df, Subtype == input$subtype)
      }
      df
    })

    #' Data for the correlation plot
    corr_two_data <- reactive({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset"),
        need(input$gene2 != "", "Please enter Gene 2")%then%
          need(input$gene2 %in% names(exprs()),"Gene not available for this dataset"),
        # Trying to avoid an error when switching datasets in case the choosen histology is not available.
        need(input$histology %in% c("All",histo()), FALSE)
      )
      data <- corr_data()[,c("Sample", "Histology", "Subtype", input$gene, input$gene2)]
      data <- data[complete.cases(data[,input$gene]),]
    })

    #' Generate the correlation plot
    output$corrPlot <- renderPlot({
      myCorggPlot(corr_two_data(), input$gene, input$gene2, colorBy = color_by(), separateBy = separate_by())
    }, width = function()input$cor_plot_width, height = function()input$cor_plot_height)

    #' Generate a summary of the correlation test
    output$corrTest <- renderTable({
      req(input$gene != input$gene2)
      df <- corr_two_data()
      if (separate_by() != "none") {
        df <- df %>% select_(group = separate_by(), gene1 = input$gene, gene2 = input$gene2)
        more_than_three <- df %>% group_by(group) %>% count(group) %>% filter(n>3)  # to drop levels with less than 3 elements
        df <- df %>% filter(group %in% more_than_three$group) %>% group_by(group)
      } else {
        df <- df %>% select_(gene1 = input$gene, gene2 = input$gene2)
      }
      cor <- df %>% do(tidy(cor.test(~ gene1 + gene2, data =., use = "complete.obs", method = tolower(input$statCorr))))
    }, digits = ifelse(input$statCorr == "Pearson",c(2,2,2,-1,2,2,2),c(2,2,2,-1)))

    #' Table with the data used for the correlation plot
    output$corrDataTable <- renderDataTable({
      data <- corr_two_data()
      data_table(data)
    })

    #' Download the corrPlot
    output$downloadcorrPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$gene, "_", input$dataset, "_corrPlot.", download_plot_file_type())
      },
      content = function(file) {
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height =  download_plot_height())
        myCorggPlot(corr_two_data(), input$gene, input$gene2, color_by(), separate_by())
        dev.off()
      }
    )

    #' Genes for multiple genes correlation
    corr_genes <- reactive({
      if(input$genelist_corr != "") {
        corr_genes <- toupper(unlist(strsplit(input$genelist_corr, split = " ")))
      } else {
        corr_genes <- input$corrGene
      }
      corr_genes
    })

    #' Multiple genes correlation
    corr_many_data <- reactive({
      validate(
        # Need two or more genes
        need(length(corr_genes()) > 1, "Please enter 2 or more genes in the panel on the left")%then%
          need(any(corr_genes() %in% names(exprs())),"Gene not available for this dataset")%then%
          # Trying to avoid an error when switching datasets in case the choosen histology is not available.
          need(input$histology %in% c("All",histo()) & input$subtype %in% c("All",subtype()), FALSE)
      )
      corr_genes <- intersect(corr_genes(), names(corr_data())) # select only the genes available for the specific dataset
      data <- corr_data()[ ,c("Sample", "Histology", "Subtype", corr_genes)]
      data <- data[complete.cases(data[,corr_genes]),]
      data
    })

    corr_many_genes <- reactive({
      corr_many_data()[,4:length(corr_many_data())]
    })

    #' Generate the pairs plot
    output$pairsPlot <- renderPlot({
      ggpairs(corr_many_genes(), mapping = aes(alpha=0.5),lower=list(continuous="smooth"))
    }, height = function(){150*length(corr_many_genes())}, width = function(){200*length(corr_many_genes())})

    #' Table with the data used for the pairs plot
    output$corrPairsDataTable <- renderDataTable({
      data <- corr_many_data()
      data_table(data)
    })

    output$missingGene <- renderText({
      missing <- setdiff(corr_genes(),names(exprs()))
      if(length(missing)>0) {
        paste("Note:", missing, "is not available for this datasets.","\n",sep = " ")
      } else {
        return()
      }
    })

    #' Download the pairs plot
    output$downloadpairsPlot <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(), "_", input$dataset, "_pairsPlot.", download_plot_file_type())
      },
      content = function(file) {
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height =  download_plot_height())
        print(ggpairs(corr_many_genes(), mapping = aes(alpha=0.5),lower=list(continuous="smooth")))
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
      row.names(corr.table) <- corr.table$Gene
      corr.table
    })

    #' Generate an HTML table view of the correlation table
    output$corrAllTable <- renderDataTable({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          # Not all genes are available for all the dataset
          need(input$gene %in% names(corr_data()),"Gene not available for this dataset")
      )
      data_table(corr_filter_data(), selection = 'single')
    })

    #' Generate a reactive value for the input$rows to set to NULL when the dataset change
    v <- reactiveValues(rows = NULL)
    observeEvent(input$corrAllTable_rows_selected, {
      v$rows <- input$corrAllTable_rows_selected
    })
    observeEvent(c(dataset_input(),input$histology,input$gene,input$cor), {
      v$rows <- NULL
    })

    #' Generate the correlation plot
    output$corrAllPlot <- renderPlot({
      validate(
        need(input$gene != "", FALSE)%then%
          need(v$rows!= "","Click on a row to see the corresponding correlation plot.")
      )
      y <- as.character(corr_filter_data()[v$rows,"Gene"]) # v$rows return only the row index, need to get the gene name
      data <- corr_data()
      data <- na.omit(corr_data()[,c(input$gene, y)])
      aes_scatter <- aes_string(x = input$gene, y = y)
      ggplot(data, mapping = aes_scatter) + theme(legend.position=c(1,1),legend.justification=c(1,1)) +
        geom_point(alpha=0.5) + geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_linedraw()
    })

    #' RPPA data analysis
    rppa_data <- reactive({
      req(input$dataset %in% c("TCGA GBM","TCGA LGG","TCGA GBMLGG"))
      req(input$gene != "")
      req(input$gene %in% names(exprs()))
      samples <- intersect(row.names(rppas()),data()[,"Sample"])
      data <- rppas()[samples,]
      mRNA <- round(data()[samples,"mRNA"],2)
      rppa_data <- list(data = data, mRNA = mRNA)
    })

    #' Create a rug plot with the mRNA expression value for the manual cutoff
    output$boxRppaRNA <- renderPlot({
      req(input$rppaCut)
      mRNA <- rppa_data()[["mRNA"]]
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
      sliderInput(inputId = "rppaCut",label = "mRNA cutoff", min = min( rppa_data()[["mRNA"]]), max = max( rppa_data()[["mRNA"]]),
                  value = median(rppa_data()[["mRNA"]]), step = 0.05, round = -2)
    })


    rrppa_data_table <- reactive({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG","TCGA GBMLGG"), "RPPA data available only for TCGA datasets")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset")%then%
          need(input$rppaCut != "", FALSE)%then%
          need(input$rppaCut > sort(rppa_data()[["mRNA"]], FALSE)[2] & input$rppaCut < sort(rppa_data()[["mRNA"]], TRUE)[2], FALSE)
      )
      rppa <- rppa_data()[["data"]]
      mRNA <- rppa_data()[["mRNA"]]
      strat <- ifelse(mRNA >= input$rppaCut, c("high"),c("low"))
      strat <- factor(strat,levels = c("low", "high"))
      d <- NULL
      results <- for (i in names(rppa)) {
        prot <- rppa[, i]
        t <- try(t.test(prot~strat),silent=TRUE)
        est <- t$estimate
        p <-  t$p.value
        x <- c(est,p)
        names(x) <- c("Avg_low","Avg_high","p.value")
        d= rbind(d,t(data.frame(x)))
      }
      row.names(d) <- names(rppa)
      d <- data.frame(Protein = row.names(d), signif(d,3))
      d$adj.p.value <- signif(p.adjust(d$p, method = "bonferroni"),3)
      d <- d[order(d$p.value),]
      d
    })

    output$rppaTable <- renderDataTable({
      data_table(rrppa_data_table(), selection = 'single')
    })

    #' Generate a reactive value for the input$rows that set to NULL when the dataset change
    rp <- reactiveValues(rppa.rows = NULL)
    observeEvent(input$rppaTable_rows_selected, {
      rp$rppa.rows <- input$rppaTable_rows_selected
    })
    observeEvent(c(dataset_input(),input$histology,input$gene,input$rppaCut), {
      rp$rppa.rows <- NULL
    })

    #' Generate the RPPA box plot
    output$rppaPlot <- renderPlot({
      protein <- as.character(rrppa_data_table()[rp$rppa.rows,"Protein"])
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG", "TCGA GBMLGG"), FALSE)%then%
          need(input$gene != "", FALSE)%then%
          need(protein!= "","Click on a row to see the corresponding plots.")%then%
          need(protein %in% names(rppa_data()[["data"]]), FALSE)
      )
      mRNA <- rppa_data()[["mRNA"]]
      strat <- ifelse(mRNA >= input$rppaCut, c("high"),c("low"))
      strat <- factor(strat,levels = c("low", "high"))
      data <- data.frame(mRNA, strat, rppa = rppa_data()[["data"]][,protein])
      r <- round(cor.test(data$mRNA,data$rppa, use = "complete.obs")$estimate,3)
      p.value <- signif(cor.test(data$mRNA,data$rppa, use = "complete.obs")$p.value,3)
      p1 <- ggplot(data, aes(x=strat, y = rppa)) + geom_boxplot(outlier.size = 0, outlier.stroke = 0) +
        geom_jitter(aes(colour = strat), position = position_jitter(width = .5), size = 2, alpha = 0.5) +
        xlab(paste(input$gene, "mRNA")) + ylab(paste(protein,"RPPA score")) +
        guides(colour=FALSE) + theme_linedraw()
      p2 <- ggplot(data, aes(x=mRNA, y = rppa)) + geom_point(aes(colour = strat), alpha=0.5) +
        geom_smooth(method = "lm", se = TRUE) + geom_rug(alpha = 0.1) + theme_linedraw() +
        xlab(sprintf("%s mRNA (log2) \n\n r = %s; p value = %s", input$gene, r, p.value)) +
        ylab(paste(protein,"RPPA score")) + theme(legend.position = "none")
      grid.arrange(p1, p2, ncol=1)
    },height = 600)

    #' Genes for mutation analysis
    mut_genes <- reactive({
      if(input$genelist_mut != "") {
        mut_genes <- toupper(unlist(strsplit(input$genelist_mut, split = " ")))
      } else {
        mut_genes <- input$mutGene
      }
      mut_genes
    })

    #' Mutation call using cgdsr
    mut_data <- reactive ({
      closeAlert(session, "mutAlertId")
      mycgds = CGDS("http://www.cbioportal.org/")
      cohort <- switch(input$dataset,
                       "TCGA GBM" = "gbm_tcga_",
                       "TCGA LGG" = "lgg_tcga_",
                       "TCGA GBMLGG" = "lgggbm_tcga_pub_"
      )
      caseList <- ifelse(cohort %in% c("gbm_tcga_","lgg_tcga_"), "sequenced","all")
      ann <- getMutationData(mycgds, caseList = paste0(cohort, caseList), geneticProfile = paste0(cohort,"mutations"), genes = mut_genes())
      ann$case_id <- gsub(pattern = "-",replacement = ".",ann$case_id)
      ann$case_id <- substr(ann$case_id,0,12)
      ann <- ann[, c(3,2,6,8,13:17)] # select only useful columns
      names(ann)[1] <- "Sample"
      if(nrow(ann) < 1) {
        createAlert(session, anchorId = "mutAlert", alertId = "mutAlertId", title = "Sorry...",
                    content = paste("No mutations identified for", input$mutGene, "in the", input$dataset, "dataset"),
                    style = "danger")
      }
      mut_df <- getProfileData(mycgds, caseList = paste0(cohort, caseList), geneticProfile = paste0(cohort,"mutations"), genes = mut_genes())
      # row.names(mut_df) <- substr(row.names(mut_df),0,12)
      mut_df <- apply(mut_df,2,as.factor)
      mut_df[mut_df == "NaN"] = ""
      mut_df[is.na(mut_df)] = ""
      mut_df[mut_df != ''] = "MUT"
      mat <- as.matrix(t(mut_df))
      if(input$add_cna){
        cna <- getProfileData(mycgds, caseList = paste0(cohort,caseList), geneticProfile = paste0(cohort,"gistic"), genes = mut_genes())
        cna <- apply(cna,2,function(x)as.character(factor(x, levels = c(-2:2), labels = c("HOMDEL", "HETLOSS", "DIPLOID", "GAIN", "AMP"))))
        cna[is.na(cna)] = ""
        cna[cna == 'DIPLOID'] = ""
        comb <- data.frame(matrix(paste(as.matrix(cna),as.matrix(mut_df),sep=";"),nrow=nrow(cna), ncol=ncol(cna),
                                  dimnames =list(row.names(mut_df),colnames(cna))))
        mat <- as.matrix(t(comb))
      }
      mut_list <- list(ann = ann, mat = mat)
    })

    output$mut <- renderDataTable({
      req(input$dataset %in% c("TCGA GBM","TCGA LGG","TCGA GBMLGG"))
      req(mut_genes() != "")
      req(any(mut_genes() %in% names(exprs())))
      mut <- mut_data()[["ann"]]
      datatable(mut, rownames = FALSE, extensions = c("FixedColumns", 'Buttons'), filter="top",
                options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE))
    })

    #' Oncoprint
    oncoprint_height <- reactive({
      if(length(mut_genes()) <= 1) {
        # if(input$column_barplot) {
        #   h <- 120
        # } else {
          h <- 65
        # }
      } else {
        h <- length(mut_genes())*55
      }
      h
    })
    output$oncoprint <- renderPlot({
      validate(
        need(input$dataset %in% c("TCGA GBM","TCGA LGG","TCGA GBMLGG"), "Mutations data available only for TCGA datasets")%then%
          need(mut_genes() != "", "Please, enter a gene or a list of genes in the panel on the left")%then%
          need(any(mut_genes() %in% names(exprs())),"Gene(s) not available for this dataset")%then%
          need(nrow(mut_data()[["ann"]])>1,FALSE)
      )
      mat <- mut_data()[["mat"]]
      alt <- apply(mat,1,function(x)strsplit(x, ";"))
      alt <- unique(unlist(alt))
      alt <- alt[which(alt !="")]
      alt <- c("background",alt)
      alt <- intersect(names(alter_fun_list), alt)
      alt_fun_list <- alter_fun_list[alt]
      col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue3","HETLOSS" = "cadetblue1", "GAIN"= "pink")
      col <- col[alt]
      oncoPrint(mat = mat, alter_fun = alt_fun_list,
                get_type = function(x) strsplit(x, ";")[[1]],
                col = col, remove_empty_columns = input$hide_cases, column_title = " ",
                top_annotation_height = unit(0, "points"),
      # top_annotation_height = unit(ifelse(input$column_barplot,55,0), "points"),          
      show_row_barplot = input$row_barplot, barplot_ignore = "NA")
    }, height = oncoprint_height)

    #' Differential Expression
    observe({
      req(input$dataset)
      pData <- rmNA(pDatas())
      colnames <- names(pData)[!sapply(pData, is.numeric)] # remove muneric categories
      colnames <- colnames[colnames %notin% c("Sample","ID","Patient_ID","Sample_ID", "Matching.sample", "Therapy_Class","title")]
      updateSelectInput(session, inputId = "colorSideHeatmap", choices = colnames)
    })

    de_data <- reactive({
      req(input$tab1 ==7)
      validate(
        need(input$dataset %notin% rnaseq_datasets, "Differential gene expression analysis currently available only for Microarray data")%then%
          need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()),"Gene not available for this dataset") # Not all genes are available for all the dataset
      )
      df <- corr_data() # data filtetered by histology and subtype
      # strat <- get_cutoff(df[,input$gene], input$DEcutoff, NULL)
      strat <- na.omit(data.frame(Sample = df$Sample, strat = get_cutoff(df[,input$gene], input$DEcutoff, NULL)))
      strat$strat <- factor(strat$strat,levels = c("low","high"))
      eset <- as.matrix(t(df[df$Sample %in% strat$Sample,-c(1:8)]))
      design <- model.matrix(~strat$strat)
      fit <- eBayes(lmFit(eset,design))
      topT <- topTable(fit, coef = 2, number = Inf, lfc = input$lfc, p.value = as.numeric(input$pvalueDE), resort.by="logFC") # genes with at least 1.4 fold change and a pvalue<.05
      if(dim(topT)[1]==0) {
        createAlert(session, anchorId = "DEAlert", alertId = "DEAlertId", title = "Sorry...",
                    content = paste("No differentially regulated genes were identified using the current settings.
                                    Consider to change either the cutoff, the Log2 fold change filter or the p value"),
                    style = "danger")
      } else {
        closeAlert(session, "DEAlertId")
      }
      esetSel <- eset[row.names(topT), ]
      if (input$pDataHeatmap) {
        req(input$colorSideHeatmap %in% names(pDatas()))
        pData <- pDatas()[pDatas()$Sample %in% strat$Sample, ]
        n <- nlevels(pData[,input$colorSideHeatmap])
        colors <- data.frame(levels(as.factor(pData[,input$colorSideHeatmap])), color = brewer_pal(palette = "Dark2")(n))
        names(colors)[1] <- paste0(input$colorSideHeatmap)
        dataColor <- merge(pData, colors, all=T)
        dataColor <- dataColor[order(dataColor$Sample),]
        heat_data <- list(eset= esetSel, topTable = topT, sideColors=dataColor, colors = colors)
      } else {
        heat_data <- list(eset= esetSel, topTable = topT)
      }
      heat_data
    })

    output$DEheatmap <- renderPlot({
      if(dim(de_data()[["topTable"]])[1]==0){
        return()
      }
      validate(
        need(dim(de_data()[["topTable"]])[1]>1, message = "Only one gene was differentially expressed with the current settings, heatmap will not be rendered"),
        errorClass = "dangerous" # errorClass does not work
      )
      if(input$pDataHeatmap) {
        req(input$colorSideHeatmap %in% names(pDatas()))
        dataColor <-de_data()[["sideColors"]]
        colors <- de_data()[["colors"]]
        col <- as.character(dataColor$color)
        heatMap <-heatmap3(de_data()[["eset"]], Rowv = NA, scale = "row", col = colorRampPalette(c("green","black", "red"))(1024), margins = c(0.5,11),
                           labCol=NA, cexRow = 1.25, ColSideColors = col, ColSideLabs = paste0(input$colorSideHeatmap))
        legend(0.85, 1, legend = colors[,input$colorSideHeatmap], fill = as.character(colors[,"color"]),
               border = FALSE, bty = "n", y.intersp = 1, cex = 1, title = paste0(input$colorSideHeatmap))
      } else {
        heatmap3(de_data()[["eset"]], Rowv = NA, scale = "row", col = colorRampPalette(c("green","black", "red"))(1024), margins = c(0,11),labCol=NA, cexRow = 1.25)
      }
    })

    output$DETable <- renderDataTable({
      validate(
        need(dim(de_data()[["topTable"]])[1]>0, message = "No differentially regulated genes were identified using the current setting.
             Consider to change either the cutoff, the Log2 fold change filter or the p value")
      )
      data <- de_data()[["topTable"]]
      data_table(data)
    })

    #' Gene Onthology
    entrez <- eventReactive(input$goGO | input$goKegg ,{
      validate(
        need(dim(de_data()[["topTable"]])[1]>5,
             "Less than 5 genes were differentially expressed with the current settings, enrichment analysis will not be performed")
      )
      entrez <-  suppressWarnings(bitr(row.names(de_data()[["topTable"]]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
      entrez <- merge(entrez,de_data()[["topTable"]], by.x="SYMBOL",by.y="row.names",all.x=T)
      entrez
    })


    go <- reactiveValues(GO = NULL, KEGG = NULL)
    observeEvent(input$goGO, {
      go$GO <- input$goGO
    })
    observeEvent(input$goKegg, {
      go$KEGG <- input$goKegg
    })
    observeEvent(c(de_data(), input$gene, input$ont, input$pvalueCutoff, input$qvalueCutoff), {
      go$GO <-  NULL
    })
    observeEvent(c(de_data(), input$gene, input$pvalueCutoffKegg, input$qvalueCutoffKegg), {
      go$KEGG <-  NULL
    })

    enrich_GO <- eventReactive(go$GO ,{
      ego <- enrichGO(gene = entrez()$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = input$ont, pAdjustMethod = "BH",
                      pvalueCutoff  = input$pvalueCutoff, qvalueCutoff  = input$qvalueCutoff, readable=TRUE)
      if(length(ego@geneInCategory)==0) {
        createAlert(session, anchorId = "goAlert", alertId = "goAlertId", title = "Sorry...",
                    content = paste("No categories identified, try to change p/q value cutoff"),
                    style = "danger")
      } else {
        closeAlert(session, "goAlertId")
      }
      ego
    })

    observe({
      l <- length(enrich_GO()@geneInCategory)
      if(l>50) l <- 50
      updateSliderInput(session,inputId = "showCategory", max = l, value = 10)
    })

    output$enrichGOPlot <- renderPlot({
      validate(
        need(!is.null(go$GO), message = "Please press the Submit button to initiate the analysis")%then%
          need(length(enrich_GO()@geneInCategory)>0, FALSE)
      )
      dotplot(enrich_GO(), showCategory= input$showCategory)
    })

    output$enrichGOMap <- renderPlot({
      req(length(enrich_GO()@geneInCategory)>0)
      # cnetplot(enrich_GO(), foldChange=entrez()$logFC,showCategory= input$showCategory)
      enrichMap(enrich_GO(), n=input$showCategory)
    })

    output$enrichGOTable <- renderDataTable({
      req(length(enrich_GO()@geneInCategory)>0)
      data_table(summary(enrich_GO()))
    })

    #' KEGG
    enrich_Kegg <- eventReactive(go$KEGG,{
      eKegg <- enrichKEGG(gene = entrez()$ENTREZID, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff  = input$pvalueCutoffKegg,
                          qvalueCutoff  = input$qvalueCutoffKegg)
      if(length(eKegg@geneInCategory)==0) {
        createAlert(session, anchorId = "keggAlert", alertId = "keggAlertId", title = "Sorry...",
                    content = paste("No categories identified, try to change p/q value cutoff"),
                    style = "danger")
      } else {
        closeAlert(session, "keggAlertId")
      }
      eKegg
    })

    observe({
      l <- length(enrich_Kegg()@geneInCategory)
      updateSliderInput(session,inputId = "showCategoryKegg", max = l, value = 10)
    })

    output$enrichKeggPlot <- renderPlot({
      validate(
        need(!is.null(go$KEGG), message = "Please press the Submit button to initiate the analysis")%then%
          need(length(enrich_Kegg()@geneInCategory)>0, FALSE)
      )
      dotplot(enrich_Kegg(), showCategory= input$showCategoryKegg)
    })

    output$enrichKeggMap <- renderPlot({
      req(length(enrich_Kegg()@geneInCategory)>0)
      enrichMap(enrich_Kegg(),n=input$showCategoryKegg)
    })

    output$enrichKeggTable <- renderDataTable({
      req(length(enrich_Kegg()@geneInCategory)>0)
      eKegg <- setReadable(enrich_Kegg(), OrgDb = "org.Hs.eg.db", keytype = "ENTREZID")
      data_table(summary(eKegg))
    })

    #' Generate reports
    output$reportPlots <- renderUI({
      validate(
        need(input$gene != "", "Please, enter a gene name in the panel on the left")%then%
          need(input$gene %in% names(exprs()), "Gene not available for this dataset")
      )
      groups <- c(plotList[[input$dataset]], plot_user_selection())
      plot_output_list <- lapply(groups, function(i) {
        plot_report <- paste("plotReport", i, sep = "")
        shinydashboard::box(height = 355, title = paste0(i), width = NULL, solidHeader = TRUE, status = "primary",
                            busy(paste("Rendering", i, "boxplot")),
                            plotOutput(plot_report, height = 300)
        )
      })
      do.call(tagList, plot_output_list)
    })

    observe({
      data <- data()
      groups <- c(plotList[[input$dataset]], plot_user_selection())
      for (i in groups) {
        local({
          my_i <- i
          plot_report <- paste("plotReport", my_i, sep = "")
          output[[plot_report]] <- renderPlot({
            validate(
              need(!all(is.na(data[ ,my_i])),"Sorry,no gene data available for this group")
            )
            data <- filter(data,!is.na(data[,my_i]))
            p <- ggplot(data, mapping=aes_string(x=my_i, y = "mRNA")) + geom_boxplot(outlier.stroke = 0, outlier.size = 0) +
              geom_jitter(position = position_jitter(width = .5), size = 2, alpha = 0.5) +
              ylab("mRNA expression (log2)") + theme(axis.title.x = element_blank()) + theme(axis.title.y=element_text(vjust=1))
            stat <- data %>%
              group_by_(my_i) %>%
              summarise(Sample_count = paste0(n()," (", round(n()*100/dim(data)[1],2), "%)" ), # prop table
                        median = round(median (mRNA, na.rm=T),2), mad = round(mad(mRNA, na.rm=T),2),
                        mean = round(mean(mRNA, na.rm=T),2), sd = round(sd(mRNA, na.rm=T),2)) %>%
              data.frame()
            row.names(stat) <- stat[,my_i]
            tot <- data %>%
              summarise(Sample_count = n(),
                        median = round(median (mRNA, na.rm=T),2), mad = round(mad(mRNA, na.rm=T),2),
                        mean = round(mean(mRNA, na.rm=T),2), sd = round(sd(mRNA, na.rm=T),2))
            stat <- stat[,-1]
            stat <- rbind(stat,TOTAL = tot)
            t <- tableGrob(stat)
            # , gp = gpar(fontsize=14),row.just = "right", core.just = "right")
            #             t <- tableGrob(stat, theme = ttheme_default(core = list(fg_params = list(hjust=1, x=0.9)),
            #                                                         rowhead = list(fg_params=list(hjust=1, x=0.95))))
            grid.arrange(p, t, ncol = 2)
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
      datatable(data, rownames = FALSE, extensions = c("FixedColumns", 'Buttons'),
                options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE))
    })

    #' Generate a graphic summary of the dataset, using googleVis
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

    # observe ({
    #   req(input$dataset)
    #   data <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status")]
    #   data <- rmNA(data)
    #   data <- subset(data, Histology != "Non-tumor")
    #   data$Histology <- droplevels(data$Histology)
    #   # Call renderChart for each one.
    #   for (i in names(data)) {
    #     local({
    #       my_i <-i
    #       plotname <- paste("plot", my_i, sep="")
    #       output[[plotname]] <- renderGvis({
    #         plotData <- data.frame(table(data[, my_i]))
    #         pie <-gvisPieChart(labelvar = "Var1", numvar = "Freq", data = plotData,
    #                            options = list(width=400, height=300, pieSliceText='label', chartArea.left = 1,
    #                                           colors = "['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628']"), chartid= plotname) # title = my_i,
    #         return(pie)
    #       })
    #     })
    #   }
    # })

    #' Generate survival groups stratified by Histology, etc.
    output$survPlots <- renderUI({
      validate(
        need(input$dataset %notin% no_surv_dataset, "Sorry, no survival data are available for this dataset")
      )
      df <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status","survival", "status")]
      df <-  rmNA(df) # Removing unavailable (all NA) groups
      groups <- names(df)[names(df) %notin% c("Sample","status","survival")]
      plot_output_list <- lapply(groups, function(i) {
        plot_surv_name <- paste("plotSurv", i, sep = "")
        plotOutput(plot_surv_name, height = 300, width = 400)
      })
      do.call(tagList, plot_output_list)
    })

    observe({
      req(input$dataset)
      df <- exprs()[ ,c("Histology", "Grade", "Recurrence", "Subtype", "CIMP_status","survival", "status")]
      df <- rmNA(df)
      if ("Histology" %in% names(df)) {
        df <- subset(df, Histology != "Non-tumor")
        df$Histology <- droplevels(df$Histology)
      } # Some of the pediatric datasets don't have histology
      groups <- names(df)[names(df) %notin% c("Sample","status","survival")]
      for (i in groups) {
        local({
          my_Survi <- i
          plot_surv_name <- paste0("plotSurv", my_Survi)
          output[[plot_surv_name]] <- renderPlot({
            req(input$dataset %notin% no_surv_dataset)
            df1 <- na.omit(data.frame(status = df[ ,"status"], time = df[ ,"survival"], strata = df[ ,my_Survi]))
            df1$strata <- droplevels(df1$strata)
            fit <- survfit(Surv(time, status == 1) ~ strata, data = df1)
            ggsurvplot(fit, legend = c(0.75,0.75), surv.scale = "percent", ylab = "% Surviving", legend.labs = levels(df1$strata), color = "red",
                       xlab = "Survival time (Months)", main = paste0("\n",my_Survi), legend.title = "", font.legend = 12, palette = "Set1")
          })
          plotname <- paste("plot", my_Survi, sep="")
          output[[plotname]] <- renderGvis({
            plotData <- data.frame(table(df[, my_Survi]))
            pie <- gvisPieChart(labelvar = "Var1", numvar = "Freq", data = plotData, chartid= plotname, options = list(width = 400, height = 300, pieSliceText = 'label', chartArea.left = 1,
                                colors = "['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#a6cee3','#006d2c','#810f7c','#02818a']")) # title = my_i,
            return(pie)
          })
        })
      }
    },priority = 100)

    #' Reactivity required to display download button after file upload
    output$finishedUploading <- reactive({
      if (is.null(input$upFile))
      { 0 } else { 1 }
    })
    outputOptions(output, 'finishedUploading', suspendWhenHidden=FALSE)

    #' Reactive function to generate SVM subtype call to pass to data table and download handler
    #' IMPORTANT #http://stackoverflow.com/questions/15503027/why-are-probabilities-and-response-in-ksvm-in-r-not-consistent
    svm_call <- eventReactive (input$goSvm | input$goSub3,{
      # input$upFile will be NULL initially. After the user selects
      # and uploads a file, it will be a data frame with 'name',
      # 'size', 'type', and 'datapath' columns. The 'datapath'
      # column will contain the local filenames where the data can
      # be found.
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,stringsAsFactors=FALSE)
      if (input$tumorType == "gbm") {
        train <- gbm.tcga[["expr"]]
        train <- subset(train,Histology=="GBM")
        Training <- train$Subtype # Expression subtype
        subtypes <- c("Classical","Mesenchymal","Proneural")
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
      set.seed(12345)
      svm <- ksvm(t(df.train), Training, cross=10, kernel="vanilladot", family="multinomial", prob.model=TRUE, scale=FALSE)
      svm.subtype.call <- as.matrix(predict(svm, t(df.learn)))
      prob <- as.matrix(predict(svm, t(df.learn), type="probabilities"))
      svm_call <- data.frame(Sample = rownames(upData), svm.subtype.call, round(prob,3))
      svm_call
    })

    #' Rerndering the subtype call as a data table
    output$svm <- renderDataTable({
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goSvm != 0,'Please press "Submit SVM"')
      )
      data_table(svm_call())
    })

    #' Reactive function to generate k-nearest neighbour subtype call to pass to data table and download handler
    knn_call <- eventReactive(input$goKnn | input$goSub3, {
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      if (input$tumorType == "gbm") {
        train <- gbm.tcga[["expr"]]
        train <- subset(train,Histology=="GBM")
        Training <- train$Subtype # Expression subtype
        subtypes <- c("Classical","Mesenchymal","Proneural")
        k <- 15
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
      set.seed(12345)
      pred <- knn3Train(train = train.exp[,genes], test = learn.exp[,genes], cl = Training, k = k, prob=TRUE)
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
      data_table(knn_call())
    })

    #' Reactive function to generate ssGSEA call to pass to data table and download handler
    gsva_call <- eventReactive(input$goGsva | input$goSub3, {
      validate(
        need(input$tumorType == "gbm","ssGSEA analysis currently available for GBM samples only"))
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE, row.names = 1)
      exprs <- data.frame(t(upData))
      # make.gct(exprs) # generate .gct
      # suppressWarnings(run.ssgsea.GBM(profile_data_file = "www/ssgsea_temp/exprs.gct", number_perms = 10))
      set.seed(12345)
      suppressWarnings(run.ssgsea.GBM(exprs, number_perms = 100))
    })

    #' Rerndering the subtype call as a data table
    output$gsva <- renderDataTable({
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goGsva != 0,'Please press "Submit ssGSEA"')
      )
      data_table(gsva_call(), rownames = T)
    })

    #' Reactive function to generate the 3 subtype calls to pass to data table and download handler
    sub3_call <- eventReactive (input$goSub3, {
      validate(
        need(input$tumorType == "gbm","'3-way' analysis currently available for GBM samples only")
      )
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,stringsAsFactors=FALSE)
      sub3 <- data_frame(Sample = upData$Sample, svm_call = svm_call()[,"svm.subtype.call"],
                         knn_call = as.character(knn_call()[,"knn.subtype.call"]), gsea_call = as.character(gsva_call()[,"gsea.subtype.call"]),
                         equal.call = ifelse(svm_call==knn_call & knn_call==gsea_call, "TRUE", "FALSE"))
      sub3$majority.call <- apply(sub3[,2:4],1,maj) # compare calls
      sub3 <- data.frame(sub3) # the dataframe generated with dplyr gives issue when sorting the html table
    })

    #' Rerndering the subtype call as a data table
    output$sub3 <- renderDataTable({
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(input$goSub3 != 0, 'Please press "Submit 3-Way"')
      )
      data_table(sub3_call())
    })

    output$sub3Plot <- renderPlot({
      data <- sub3_call()[,1:4]
      data$Sample <- factor(data$Sample, levels=(data$Sample)[order(data$svm_call)]) # sort the sample on the svm_call
      data <- melt(data, id.var = "Sample")
      names(data) <- c("Sample","Call","Subtype")
      ggplot(data, aes(Sample, Call)) + geom_tile(aes(fill = Subtype), colour = "white") + ylab("") + xlab("") +
        scale_x_discrete(expand = c(0, 0)) + theme_minimal() + scale_fill_brewer(palette = "Set1") + # scale_fill_manual(values = terrain.colors(4)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + theme(legend.position="bottom")
    }, height = 200)

    #' Reactive function to generate Estimate call to pass to data table
    est_call <- eventReactive (input$goEst,{
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
      est_call()
    }, selection = 'single', extensions = 'Buttons',
    options = list(orderClasses = TRUE, columnDefs = list(list(visible = FALSE, targets = 0)), autoWidth = TRUE)
    )

    #' Generate the Purity plot
    output$purityPlot <- renderPlot({
      validate(
        need(est_call(), FALSE)%then%
          need(input$platformEst == "affymetrix", "Sorry, the purity plot for Agilent and RNAseq have not yet been impemented") %then%
          need(input$estScore_rows_selected != "","Click on a row to see the corresponding purity plot.")
      )
      plotPurity(est_call(), est_call()[input$estScore_rows_selected,"Sample"], platform = input$platformEst)
    })

    #' Gene sets in msigdb.v5
    updateSelectizeInput(session, inputId = "genesets_msigdb", choices = row.names(msigdb.v5), server = TRUE)

    goDec <- reactiveValues(Submit = NULL)
    observeEvent(input$goDec, {
      goDec$Submit <- input$goDec
    })
    observeEvent(c(input$upFile, input$genesets_msigdb, input$geneListDec,input$min.sz,input$max.sz), {
      goDec$Submit <-  NULL
    })

    #' Reactive function to generate Deconvolute scores to pass to data table
    deconv_call <- eventReactive(goDec$Submit,{
      inFile <- input$upFile
      upData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, stringsAsFactors=FALSE)
      row.names(upData) <- upData[,"Sample"]
      exprs <- data.frame(t(upData[,-1]),check.names=FALSE)
      platformDec <- ifelse(input$platformDec == "rnaseq",TRUE, FALSE)
      if(input$geneListDec == "MSigDB gene sets"){
        gene_list <- msigdb.v5[input$genesets_msigdb,-1]
        gene_list <- as.list(data.frame(t(gene_list),stringsAsFactors = FALSE))
        # gene_list <- getBroadSets(asBroadUri(input$genesets_msigdb))
      } else {
        gene_list <- switch(input$geneListDec,
                            "Newman et al. 2015" = LM22_gene_set_list,
                            "Engler et al. 2012" = engler_gene_set_list,
                            "Bindea et al. 2013" = galon_gene_set_list)
      }
      set.seed(1234)
      gsva_results <- gsva(expr=as.matrix(exprs), gset.idx.list = gene_list, method="ssgsea", rnaseq = platformDec, parallel.sz = 0,
                           min.sz = input$min.sz, max.sz= input$max.sz, verbose= FALSE)
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
      if(input$geneListDec == "MSigDB gene sets"){
        validate(
          need(input$genesets_msigdb != "", "Please specify MSigDB gene set(s)")%then%
            need(length(input$genesets_msigdb)>1, "Two or more gene sets are required to render the heatmap")
        )
      }
      validate(
        need(!is.null(input$upFile),"Please upload the dataset to be analyzed")%then%
          need(!is.null(goDec$Submit),'Please press "Submit Deconvolute"')
      )
      data <- deconv_call()[["results"]]
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
          need(!is.null(goDec$Submit),'Please press "Submit Deconvolute"')
      )
      data <- deconv_call()[["results"]]
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
          need(!is.null(goDec$Submit),'Please press "Submit Deconvolute"')
      )
      datatable(deconv_call()[["scores"]], rownames = FALSE, extensions = c("FixedColumns", 'Buttons'),
                options = list(scrollX = TRUE, scrollCollapse = TRUE, orderClasses = TRUE, autoWidth = TRUE)
      )
    })

    #' Download expression data
    output$downloadExpressionData <- downloadHandler(
      filename = function() {
        paste(Sys.Date(), input$dataset, "expression.txt", sep = "_")
      },
      content = function(file) {
        write.table(exprs()[,-c(2:8)], file, sep = "\t", row.names = FALSE)
      }
    )

    #' Download pData
    output$downloadpData <- downloadHandler(
      filename = function() {
        paste(Sys.Date(), input$dataset, "pheno.txt", sep = "_")
      },
      content = function(file) {
        write.table(pDatas(), file, sep = "\t", row.names = FALSE)
      }
    )
  })
