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
                 set <- read.FCS(fcsFiles$datapath) })
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
        return(timestep)
    })


    ## order fcs expression according acquisition time
    ordFCS <- reactive({
        if(is.null(set()))
            return(NULL)
        ordFCS <- ord_fcs_time(set(), timeChannel())
        return(ordFCS)
    })


    ## flow rate time slider UI
    output$timeSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        flowRateData <- cellCheck()[[1]]
        mint <- min(flowRateData$frequencies[,3]) - 0.1
        maxt <- max(flowRateData$frequencies[,3]) + 0.1
        sliderInput("timeSlider", strong("Time cut:"),
                    min = mint, max = maxt, value = c(mint, maxt), step = 0.1)
    })


    output$rateSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        flowRateData <- cellCheck()[[1]]
        minc <- min(flowRateData$frequencies[,4]) - 10
        maxc <- max(flowRateData$frequencies[,4]) + 10
        sliderInput("rateSlider", strong("Flow rate cut:"),
                    min = minc, max = maxc, value = c(minc, maxc), step = 0.1)
    })

    output$signalBinSlider <- renderUI({
        if(is.null(set()) || is.null(cellCheck()))
            return(NULL)
        flowSignalData <- cellCheck()[[2]]
        maxb <- nrow(flowSignalData$exprsBin)
        sliderInput("signalBinSlider", strong("Flow rate cut:"), width = "90%",
                    min = 1, max = maxb, value = c(1, maxb), step = 1)
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
        numericInput("signalBinSize", label = h5("Signal bin size (event number)"),
                     value = optSize, min = 1, max = maxSize)
    })


    ## cell quality check
    cellCheck <- reactive({
        if(is.null(ordFCS()))
            return(NULL)

        flowRateData <- flow_rate_bin(ordFCS(), second_fraction = input$timeLenth,
                                      timeCh = timeChannel(), timestep = timeStep())

        flowSignalData <- flow_signal_bin(ordFCS(), channels = NULL, binSize = input$signalBinSize,
                                          timeCh = timeChannel(), timestep = timeStep())

        flowMarginData <- flow_margin_check(ordFCS())

        res <- list(flowRateData, flowSignalData, flowMarginData)
        return(res)
    })



    ## plot
    output$flowRatePlot <- renderPlot({
        if(is.null(ordFCS()) || is.null(cellCheck()))
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

        FlowRateQC <- flow_rate_check(cellCheck()[[1]], input$rateSlider[1], input$rateSlider[2],
                                      input$timeSlider[1], input$timeSlider[2])
        FlowSignalQC <- flow_signal_check(cellCheck()[[2]], input$signalBinSlider[1], input$signalBinSlider[2])

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

        good_sub_exprs <- sub_exprs[goodCellIDs, ]
        goodfcs <- flowFrame(exprs = good_sub_exprs, parameters = params, description = keyval)

        bad_sub_exprs <- sub_exprs[badCellIDs, ]
        badfcs <- flowFrame(exprs = bad_sub_exprs, parameters = params,description = keyval)

        return(list(totalCellNum, totalBadPerc, goodfcs, badfcs,
                    flowRatePerc, flowSignalPerc, flowMarginPerc))
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
        paste0("Percentage of low-Q events in flow rate check: ", round(checkRes()[[5]]*100,2), "%")
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
      file_ext <- description(ordFCS())$FIL
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

})


