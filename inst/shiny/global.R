# FlowQ test
require(RColorBrewer)
require(flowCore)
require(reshape2)
require(plyr)
require(scales)
require(ggplot2)

# Guess which channel captures time in a exprs, flowFrame or flowset
findTimeChannel <- function(xx, strict = FALSE) {
    time <- grep("^Time$", colnames(xx), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(time)) {
        if (is(xx, "flowSet") || is(xx, "ncdfFlowList"))
            xx <- exprs(xx[[1]]) else if (is(xx, "flowFrame"))
                xx <- exprs(xx)
            cont <- apply(xx, 2, function(y) all(sign(diff(y)) >= 0))
            time <- names(which(cont))
    }
    if (!length(time) && strict)
        stop("Unable to identify time domain recording for this data.\n", "Please define manually.",
             call. = FALSE)
    if (length(time) > 1)
        time <- character(0)
    return(time)
}

# Check if the Fcs file is ordered according to time otherwise it order it.
ord_fcs_time <- function(x, timeCh= "Time"){
    xord <- order(exprs(x)[, timeCh])

    if( !identical(xord, 1:nrow(x)) ){
        warning(paste0("Expression data in the file ", basename(keyword(x)$FILENAME),
                       " were not originally ordered by time."))
        params <- parameters(x)
        keyval <- keyword(x)
        sub_exprs <- exprs(x)[xord, ]
        newx <- flowFrame(exprs = sub_exprs, parameters = params,
                          description = keyval)
        return(newx)
    }else{
        return(x)
    }
}



flow_rate_bin <- function(x, second_fraction = 0.1, timeCh = "Time", timestep = 0.34){

    xx <- exprs(x)[, timeCh]
    idx <- c(1:nrow(x))
    lenx <- length(xx)                                                  # num of time ticks

    endsec <- ceiling(timestep * max(xx))                 # total seconds of the experiment
    tbins <- seq(0, endsec/timestep, by = second_fraction/timestep)             # time bins
    secbin <- seq(0, endsec, by = second_fraction)               # bin expressed in seconds
    minbin <- round(secbin/60, 3)                                # bin expressed in minutes
    tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)  # number of events per time bin

    nrBins <- length(tbins) - 1
    expEv <- lenx/(nrBins)           # median(tbCounts) # expected number of events per bin
    binID <- do.call(c, mapply(rep, x = 1:length(tbCounts), times = tbCounts,
                               SIMPLIFY = FALSE))

    if (length(idx) != length(binID))
        stop("length of cell ID not equal length of bin ID")

    flowRateData <- list(frequencies = cbind(tbins, minbin, secbin, tbCounts),
                         cellBinID = data.frame(cellID = idx, binID = binID),
                         info = data.frame(second_fraction = second_fraction,
                                           expFrequency = expEv, bins = nrBins))
    return(flowRateData)
}


flow_rate_plot <- function(flowRateData, lowerRateThres, upperRateThres,
                           lowerTimeCut, UpperTimeCut) {

    frequencies <- as.data.frame(flowRateData$frequencies)
    second_fraction <- flowRateData$info$second_fraction
    short_period <- quantile(frequencies$secbin, seq(0,1, length.out = 4))
    long_period <- quantile(frequencies$secbin, seq(0,1, length.out = 3*10 + 1))

    ## flow rate graph(frg)
    frg <- ggplot(frequencies, aes_string(x="secbin", y="tbCounts")) + geom_line(colour="red") +
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           text=element_text(size = 14))
    frg <- frg + geom_vline(xintercept=short_period, color="gray60") +
        geom_vline(xintercept=long_period, linetype = 2, color="gray60")

    frg <- frg + geom_hline(yintercept=c(lowerRateThres, upperRateThres), color="blue",
                            linetype = "longdash", size = 1.2, show_guide = TRUE)
    frg <- frg + geom_vline(xintercept=c(lowerTimeCut, UpperTimeCut), color="blue",
                            linetype = "longdash", size = 1.2, show_guide = TRUE)

    frg <- frg + labs(x= "Seconds", y= paste0("Number of cells per 1/", 1/second_fraction, " second"),
                      title= "Flow Rate Plot")

    return(frg)
}

flow_rate_check <- function(flowRateData, lowerRateThres, upperRateThres,
                           lowerTimeCut, UpperTimeCut) {

    frequencies <- as.data.frame(flowRateData$frequencies)
    cellBinID <- flowRateData$cellBinID

    goodcell_x <- which(frequencies$secbin < UpperTimeCut & frequencies$secbin > lowerTimeCut)
    goodcell_y <- which(frequencies$tbCounts < upperRateThres & frequencies$tbCounts > lowerRateThres)

    flowRateQC <- cellBinID$cellID[ cellBinID$binID %in% intersect(goodcell_x, goodcell_y) ]
    cat("flow rate high Q cells: ", length(flowRateQC), "\n")
    return(flowRateQC)
}



flow_signal_bin <- function(x, channels = NULL, binSize=500, timeCh="Time", timestep=0.34) {

    if (is.null(channels) || missing(channels) || is.na(channels)) {
        parms <- setdiff(colnames(x), timeCh)
    } else {
        if (!all(channels %in% colnames(x)))
            stop("Invalid channel(s)")
        parms <- channels
    }

    if (missing(binSize) || is.null(binSize) || is.na(binSize))
        binSize <- 500

    ### Retriving time and expression info
    exp <- exprs(x)
    timex <- exp[, timeCh]
    yy <- exp[, parms]  # channels data
    idx <- c(1:nrow(x))
    seconds <- timex * timestep
    lenSec <- length(seconds)  # num of time ticks
    uniSeconds <- unique(seconds)  # num of unique time tick
    nrBins <- floor(lenSec/binSize)  # num of bins

    if (length(uniSeconds) < nrBins || lenSec < binSize)
        stop("Improper bin size")

    cf <- c(rep(1:nrBins, each = binSize), rep(nrBins + 1, lenSec - nrBins * binSize))  # id bins
    stopifnot(length(cf) == lenSec)
    tmpx <- split(seconds, cf)
    xx2 <- sapply(tmpx, mean)       # mean of each time bin  (x axis)
    yy2 <- as.matrix(ddply(as.data.frame(yy), .(cf), colwise(median)))[, -1]

    return(list(exprsBin = cbind(timeSec = xx2, yy2), cellBinID = data.frame(cellID = idx, binID = cf),
                bins = length(unique(cf)), binSize = binSize))
}


flow_signal_plot <- function(flowSignalData, lowerBinThres, upperBinThres) {

    exprsBin <- flowSignalData$exprsBin

    binID <- 1:nrow(exprsBin)
    teCh <- grep("Time|time|Event|event", colnames(exprsBin), value = T)
    parms <- setdiff(colnames(exprsBin), teCh)
    dataORIG <- exprsBin[, parms]     # first channel is time
    data <- as.data.frame(dataORIG)
    data$binID <- binID

    longdata <- melt(data, id.vars = "binID", variable.name = "marker", value.name = "value")
    FS_graph <- ggplot(longdata, aes(x = binID, y = value, col = marker), environment = environment()) +
        geom_line() + facet_grid(marker ~ ., scales = "free") +
        labs(x = "segment Id", y = "Median Intensity value") + theme_bw() +
        theme(strip.text.y = element_text(angle = 0, hjust = 1), axis.text = element_text(size = 10),
              axis.title = element_text(size = 15), legend.position = "none") +
        scale_x_continuous(breaks= pretty_breaks(n = 10)) +
        scale_y_continuous(breaks= pretty_breaks(n = 3)) +
        geom_rect(aes(xmin = lowerBinThres, xmax = upperBinThres, ymin = -Inf, ymax = Inf), fill = "orange", linetype = 0, alpha = 0.005)

    return(FS_graph)
}


flow_signal_check <- function(flowSignalData, lowerBinThres, upperBinThres) {

    exprsBin <- flowSignalData$exprsBin
    cellBinID <- flowSignalData$cellBinID

    goodBins <- cellBinID$binID < upperBinThres & cellBinID$binID > lowerBinThres
    FlowSignalQC <- cellBinID$cellID[goodBins]

    cat("flow signal check: ", length(FlowSignalQC), "\n")
    return(FlowSignalQC)
}


flow_margin_check <- function(x,  margin_channels = NULL, side = "both") {

    if (is.null(margin_channels)) {
        teCh <- grep("Time|time|Event|event", colnames(x), value = T)
        parms <- setdiff(colnames(x), teCh)
    } else {
        if (!all(margin_channels %in% colnames(x)))
            stop("Invalid channel(s)")
        parms <- margin_channels
    }
    scatter_parms <- grep("FSC|SSC", parms, value = T)

    xx <- c(1:nrow(x))
    yy <- x@exprs[, parms]
    range <- range(x)
    lenx <- length(xx)

    ## lower check
    if (side == "lower" || side == "both") {

        out_neg_range <- apply(yy, 2, function(x) {
            neg <- which(x < 0)
            # Zscores <- (0.6745*(x[neg] + median(x[neg])))/mad(x[neg]) ## it
            # calculates the Zscore outneg <- neg[which(Zscores < -3.5)]
            min_value <- (-3.5 * mad(x[neg]) + (0.6745 * median(x[neg])))/0.6745  # -3.5 is the default threshold
            if (is.na(min_value)) {
                min(x) - 1
            } else {
                min_value
            }
        })
    }

    # n. bad cells for each channel
    if (side == "lower" || side == "both") {
        neg_bad_len <- sapply(parms, function(x) length(xx[yy[, x] <= out_neg_range[x]]))
    }
    if (side == "upper" || side == "both") {
        pos_bad_len <- sapply(parms, function(x) length(xx[yy[, x] >= range[2,
                                                                            x]]))
    }

    # badcellIDs
    if (side == "lower" || side == "both") {
        lowID <- do.call(c, lapply(parms, function(ch) {
            xx[yy[, ch] <= out_neg_range[ch]]
        }))
        if(length(scatter_parms) != 0){   ### check for values less than 0 in scatter parameters
            minSc <- apply(yy[,scatter_parms], 1, function(x){
                min(x)
            })
            low_scatter_ID <- which(minSc < 0)
            lowID <- unique(c(lowID, low_scatter_ID))
        }
    }
    if (side == "upper" || side == "both") {
        upID <- do.call(c, lapply(parms, function(ch) {
            xx[yy[, ch] >= range[2, ch]]
        }))
    }

    if (side == "lower") {
        summary_bad_cells <- data.frame(lower_range = c(neg_bad_len,
                                                        total_SUM = length(lowID), total_UNIQUE = length(unique(lowID))))
        bad_lowerIDs <- unique(lowID)
        bad_upperIDs <- NULL
        badCellIDs <- unique(lowID)
    } else if (side == "upper") {
        summary_bad_cells <- data.frame(upper_range = c(pos_bad_len,
                                                        total_SUM = length(upID), total_UNIQUE = length(unique(upID))))
        bad_lowerIDs <- NULL
        bad_upperIDs <- unique(upID)
        badCellIDs <- unique(upID)
    } else {
        summary_bad_cells <- data.frame(lower_range = c(neg_bad_len,
                                                        total_SUM = length(lowID), total_UNIQUE = length(unique(lowID))),
                                        upper_range = c(pos_bad_len,
                                                        total_SUM = length(upID), total_UNIQUE = length(unique(upID))))
        bad_lowerIDs <- unique(lowID)
        bad_upperIDs <- unique(upID)
        badCellIDs <- unique(c(lowID,upID))
    }

    goodCellIDs <- setdiff(xx, badCellIDs)

    cat("margin check:", length(goodCellIDs), "\n")

    return(list(goodCellIDs = goodCellIDs, bad_lowerIDs = bad_lowerIDs,
                bad_upperIDs = bad_upperIDs, events = lenx))
}


###  graph showing where the anomalies mostly happened
flow_margin_plot <- function(FlowMarginData, binSize = 500) {

    tot_events <- FlowMarginData$events
    bad_lowerIDs <- FlowMarginData$bad_lowerIDs
    bad_upperIDs <- FlowMarginData$bad_upperIDs

    if (missing(binSize) || is.null(binSize) || is.na(binSize))
        binSize <- 500
    nrBins <- floor(tot_events/binSize)

    cf <- c(rep(1:nrBins, each = binSize), rep(nrBins + 1, tot_events - nrBins * binSize))
    tmpx <- split(1:tot_events, cf)

    if(length(bad_lowerIDs) != 0 & length(bad_upperIDs) != 0){
        lowline <- sapply(tmpx, function(x){
            length(which(bad_lowerIDs %in% x))
        })
        upline <- sapply(tmpx, function(x){
            length(which(bad_upperIDs %in% x))
        })
        ymax <- max(lowline, upline)
        plot(lowline, type ="l", col = "blue", bty ="n",
            ylim = c(0, ymax), xlab = "segment ID",
            ylab = "Number of cells removed" )
        lines(upline, col = "red")
        legend("top", c("Negative Outliers", "Upper Margine Events"), lty = 1,bty = "n", cex = 0.7,
            col = c("blue", "red"))
    }else if( length(bad_lowerIDs) != 0 & length(bad_upperIDs) == 0){
        lowline <- sapply(tmpx, function(x){
            length(which(bad_lowerIDs %in% x))
        })
        plot(lowline, type ="l", col = "blue", bty ="n", xlab = "segment ID",
            ylab = "Number of cells removed" )
        legend("top", c("Negative Outliers"), lty = 1,bty = "n", cex = 0.7,
            col = "blue")
    }else if( length(bad_lowerIDs) == 0 & length(bad_upperIDs) != 0){
        upline <- sapply(tmpx, function(x){
            length(which(bad_upperIDs %in% x))
        })
        plot(upline, type ="l", col = "red", bty ="n", xlab = "segment ID",
            ylab = "Number of cells removed" )
        legend("top", c("Upper Margine Events"), lty = 1,bty = "n", cex = 0.7,
            col = "red")
    }
}
