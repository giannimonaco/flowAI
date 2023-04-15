## max data size
options(shiny.maxRequestSize=1024^3)

shinyServer(function(input, output, session) {

    ##-------------------------------------------------------------------------

    ## load flowset data
    set <- reactive({
        if (input$goButton == 0)
            return()
        isolate({fcsFiles <- input$fcsFiles
                 if (is.null(fcsFiles))
                     return(NULL)
                 set <- read.FCS(fcsFiles$datapath)
                 set@description$FILENAME <- fcsFiles$name})
        return(set)
    })

    ## time channel name
    timeChannel <- reactive({
        if(is.null(set()))
            return(NULL)
        x <- set()
        time <- findTimeChannel(x)
        return(time)
    })

    ## time step
    timeStep <- reactive({
        if(is.null(set()))
            return(NULL)
        word <- which(grepl("TIMESTEP", names(set()@description),
                            ignore.case = TRUE))
        timestep <- as.numeric(set()@description[[word[1]]])
        if( !length(timestep) ){
            warning("The timestep keyword was not found in the FCS file and it was set to 0.01. Graphs labels indicating time might not be correct", call. =FALSE)
            timestep <- 0.01
        }
        return(timestep)
    })


    TimeChCheck <- reactive({
        if (!is.null(timeChannel())) {
            if (length(unique(exprs(set())[, timeChannel()])) == 1){
                TimeChCheck <- "single_value"
            }else{
                TimeChCheck <- NULL 
            }
        }else{
            TimeChCheck <- "NoTime"
        }
        return(TimeChCheck)
    })

    
    ## order fcs expression according acquisition time
    ordFCS <- reactive({
        if(is.null(set()))
            return(NULL)
        if(is.null(TimeChCheck())){
          ordFCS <- ord_fcs_time(set(), timeChannel())
        }else{
          ordFCS <- set()
        }
      return(ordFCS)
    })
    

    ## signal bin size UI
    output$signalBinSize <- renderUI({
        if(is.null(set())){
            optSize <- NULL
            maxSize <- Inf
        }else{
            maxSize <- nrow(ordFCS())
            optSize <- min(max(1, floor(maxSize/100)), 500)
        }
        numericInput("signalBinSize", label = h5("Number of events per bin:"),
                     value = optSize, min = 1, max = maxSize)
    })


    ## cell quality check
    cellCheck <- reactive({
        if(is.null(ordFCS()))
            return(NULL)
        if(is.null(TimeChCheck())){
         flowRateData <- flow_rate_bin(ordFCS(), second_fraction = input$timeLenth,
                                      timeCh = timeChannel(), timestep = timeStep())
        }else{
            flowRateData <- list()
        }
        flowSignalData <- flow_signal_bin(ordFCS(), channels = NULL, 
                                          binSize = input$signalBinSize, timeCh = timeChannel(), 
                                          timestep = timeStep(), TimeChCheck = TimeChCheck() )

        flowMarginData <- flow_margin_check(ordFCS())

        res <- list(flowRateData, flowSignalData, flowMarginData)
        return(res)
    })

    
    ## flow rate time slider UI and check sliders. if they are null, a default value is returned for the QC
    sliders <- reactive({
        flowRateData <- cellCheck()[[1]]
        flowSignalData <- cellCheck()[[2]]
        return(c(
            min(flowRateData$frequencies[,3]) - 0.1,
            max(flowRateData$frequencies[,3]) + 0.1,
            min(flowRateData$frequencies[,4]) - 10,
            max(flowRateData$frequencies[,4]) + 10,
            0,
            nrow(flowSignalData$exprsBin) + 1)
            )
    })
    
    output$timeSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
            return(NULL)
        sliderInput("timeSlider", strong("Time cut:"),
            min = sliders()[1], max = sliders()[2], 
            value = c(sliders()[1], sliders()[2]), step = 0.1)
    })
    timeSlider <- reactive({
        if(is.null(input$timeSlider)){
            return(c(sliders()[1], sliders()[2]))
        }else{
            return(c(input$timeSlider[1],  input$timeSlider[2]))
        }
        
    })
    
    output$rateSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
            return(NULL)
        sliderInput("rateSlider", strong("Flow rate cut:"),
            min = sliders()[3], max = sliders()[4], 
            value = c(sliders()[3], sliders()[4]), step = 0.1)
    })
    rateSlider <- reactive({
        if(is.null(input$rateSlider)){
            flowRateData <- cellCheck()[[1]]
            return(c(sliders()[3], sliders()[4]))
        }else{
            return(c(input$rateSlider[1],  input$rateSlider[2]))
        }
        
    })
    
    output$signalBinSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        sliderInput("signalBinSlider", strong("Signal acquisition cut:"), width = "90%",
            min = sliders()[5], max = sliders()[6], 
            value = c(sliders()[5], sliders()[6]), step = 1)
    })
    signalSlider <- reactive({
        if(is.null(input$signalBinSlider)){
            return(c(sliders()[5], sliders()[6]))
        }else{
            return(c(input$signalBinSlider[1],  input$signalBinSlider[2]))
        }
    }) 
    

    ## plot
    output$flowRatePlot <- renderPlot({
        if(is.null(ordFCS()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
            return(NULL)
          flowRateData <- cellCheck()[[1]]
          frp <- flow_rate_plot(flowRateData, input$rateSlider[1], input$rateSlider[2],
                             input$timeSlider[1], input$timeSlider[2])
          print(frp)
    })

    output$flowSignalPlot <- renderPlot({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        flowSignalData <- cellCheck()[[2]]
        fsp <- flow_signal_plot(flowSignalData, input$signalBinSlider[1], input$signalBinSlider[2])
        print(fsp)
    })

    output$flowMarginPlot <- renderPlot({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        flowMarginData <- cellCheck()[[3]]
        fmp <- flow_margin_plot(flowMarginData, input$signalBinSize)
        print(fmp)
    })


    
    ## check results
    checkRes <- reactive({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)

        ordFCS <- ordFCS()
        totalCellNum <- nrow(ordFCS)
        origin_cellIDs <- 1:totalCellNum
        if(is.null(TimeChCheck())){
          FlowRateQC <- flow_rate_check(cellCheck()[[1]], rateSlider()[1], rateSlider()[2],
              timeSlider()[1], timeSlider()[2])
        }else{
          FlowRateQC <- origin_cellIDs
        }
        FlowSignalQC <- flow_signal_check(cellCheck()[[2]], signalSlider()[1], signalSlider()[2])

        if(input$checkbox[1] == TRUE){
            FlowMarginQC <- cellCheck()[[3]]$goodCellIDs
        }else{
            FlowMarginQC <- origin_cellIDs
        }

        goodCellIDs <- intersect(FlowRateQC, intersect(FlowSignalQC, FlowMarginQC))
        badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)

        flowRatePerc <- 1 - length(FlowRateQC)/length(origin_cellIDs)
        flowSignalPerc <- 1 - length(FlowSignalQC)/length(origin_cellIDs)
        flowMarginPerc <- 1 - length(FlowMarginQC)/length(origin_cellIDs)
        totalBadPerc <- length(badCellIDs)/length(origin_cellIDs)

        params <- parameters(ordFCS)
        keyval <- keyword(ordFCS)
        sub_exprs <- exprs(ordFCS)

        check_goodcellsID(goodCellIDs)
        good_sub_exprs <- sub_exprs[goodCellIDs, , drop = FALSE]
        goodfcs <- flowFrame(exprs = good_sub_exprs, parameters = params, description = keyval)

        check_badCellIDs(badCellIDs)
        bad_sub_exprs <- sub_exprs[badCellIDs, , drop = FALSE]
        badfcs <- flowFrame(exprs = bad_sub_exprs, parameters = params,description = keyval)

        tempQCvector <- cellCheck()[[2]]
        QCvector <- tempQCvector$cellBinID[,"binID"]
        QCvector[badCellIDs] <- runif(length(badCellIDs), min=10000, max=20000) 
        QCfcs <- addQC(QCvector, sub_exprs, params, keyval)
        
        return(list(totalCellNum, totalBadPerc, goodfcs, badfcs,
                    flowRatePerc, flowSignalPerc, flowMarginPerc, QCfcs))
    })

    ## summary text
    output$summaryText1 <- renderText({
        if(is.null(checkRes()))
            return(NULL)
        paste0("Total number of events: ", checkRes()[[1]])
    })

    output$summaryText2 <- renderText({
        if(is.null(checkRes()))
            return(NULL)
        paste0("Percentage of low-Q events: ", round(checkRes()[[2]]*100,2), "%")
    })

    output$flowRateSummary <- renderText({
        if(is.null(checkRes()))
            return(NULL)
        if(is.null(TimeChCheck())){
          paste0("Percentage of low-Q events in flow rate check: ", round(checkRes()[[5]]*100,2), "%")
        }else if(!is.null(TimeChCheck()) && TimeChCheck() == "NoTime"){
            "It is not possible to recreate the flow rate because the time channel is missing."
        }else if(!is.null(TimeChCheck()) && TimeChCheck() == "single_value"){
          "It is not possible to recreate the flow rate because the time channel contains a single value."
        }
    })

    output$flowSignalSummary <- renderText({
        if(is.null(checkRes()))
            return(NULL)
        paste0("Percentage of low-Q events in signal acquisition check: ", round(checkRes()[[6]]*100,2), "%")
    })

    output$flowMarginSummary <- renderText({
        if(is.null(checkRes()))
            return(NULL)
        paste0("Percentage of low-Q events in dynamic range check: ", round(checkRes()[[7]]*100,2), "%")
    })

    file_base <- reactive({
      file_ext <- description(ordFCS())$FILENAME
      file_base <- sub("^([^.]*).*", "\\1", file_ext)
      return(file_base)
    })

    ## download processed FCS files
    output$downloadGoodFCS <- downloadHandler(
        filename = function(){
          paste0(file_base(), "_HighQ.fcs")
          },

        content = function(file){
            data <- checkRes()[[3]]
            if(is.null(data)){
                return(NULL)
            }
            write.FCS(data, file)
            #tar(tarfile = file, files = tempdir)
        }
    )

    output$downloadBadFCS <- downloadHandler(
        filename = function(){
          paste0(file_base(), "_LowQ.fcs")
        },

        content = function(file){
            data <- checkRes()[[4]]
            if(is.null(data)){
                return(NULL)
            }
            write.FCS(data, file)
            #tar(tarfile = file, files = tempdir)
        }
    )

    ## download processed FCS files
    output$downloadQCFCS <- downloadHandler(
        filename = function(){
            paste(file_base(), "_QC.fcs", sep='')
        },
        
        content = function(file){
            data <- checkRes()[[8]]
            if(is.null(data)){
                return(NULL)
            }
            write.FCS(data, file)
            #tar(tarfile = file, files = tempdir)
        }
    )
    
})

