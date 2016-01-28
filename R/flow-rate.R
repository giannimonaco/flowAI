# A flowFrame object is splitted in bins that are indicate the flow
# rate over time.
# @param second_fraction # change the fraction of seconds according
# to your desire
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, values, and
# expected values.}
#
flow_rate_bin <- function(x, second_fraction = 0.1, timeCh = timeCh,
  timestep = timestep){

  xx <- exprs(x)[, timeCh]
  idx <- c(1:nrow(x))

  endsec <- ceiling(timestep * max(xx))  # total seconds of the experiment

  lenx <- length(xx)  # num of time ticks

  tbins <- seq(0, endsec/timestep, by = second_fraction/timestep)  # time bins
  secbin <- seq(0, endsec, by = second_fraction)  # bin expressed in seconds
  minbin <- round(secbin/60, 3)  # bin expressed in minutes
  nrBins <- length(tbins) - 1
  tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)  # number of events per time bin
  expEv <- lenx/(nrBins)  ##median(tbCounts) # expected number of events per bin
  binID <- do.call(c, mapply(rep, x = 1:length(tbCounts), times = tbCounts,
    SIMPLIFY = FALSE))

  if (length(idx) != length(binID))
    stop("length of cell ID not equal length of bin ID")

  timeFlowData <- list(frequencies = cbind(tbins, minbin, secbin, tbCounts),
                       cellBinID = data.frame(cellID = idx, binID = binID),
                       info = data.frame(second_fraction = second_fraction,
                       expFrequency = expEv, bins = nrBins))
  return(timeFlowData)
}


# Detection of anomalies in the flow rate using the algorithm
# implemented in the package AnomalyDetection.
flow_rate_check <- function(x, FlowRateData,
        long_period = long_period, period = period , alpha = alpha) {
  fr_frequences <- FlowRateData$frequencies
  fr_cellBinID <- FlowRateData$cellBinID
  second_fraction <- FlowRateData$info["second_fraction"]

    long_period_obs <- round(nrow(fr_frequences)/long_period)
    period_obs <- ceiling(long_period_obs/period)

  if (length(unique(fr_frequences[, 2])) == 1) {
    fr_autoqc <- NULL
  } else {
    fr_autoqc <- anomaly_detection(fr_frequences[, "tbCounts"],
                 max_anoms = 0.4, period = period_obs, direction = "both",
                 longterm_period = long_period_obs, alpha = alpha)
  }

  if (is.null(fr_autoqc) || is.null(fr_autoqc$anoms)) {
    cat("NO abnormal cell detected in time flow check!", fill = TRUE)
    badPerc <- 0
    newx <- x
    goodCellIDs <- fr_cellBinID$cellID
    badCellIDs <- NULL
  } else {
    goodCellIDs <- fr_cellBinID$cellID[!(fr_cellBinID$binID %in% fr_autoqc$anoms$index)]
    badCellIDs <- setdiff(fr_cellBinID$cellID, goodCellIDs)
    badPerc <- round(1 - (length(goodCellIDs)/nrow(fr_cellBinID)), 2)
    if (badPerc > 0.5) {
      warning(paste0(100 * badPerc, "%", " of cells are detected as abnormal cells in the flow rate check. \n"))
    } else {
      cat(paste0(100 * badPerc, "%", " of cells are detected as abnormal cells in the flow rate check. \n"))
    }

    params <- parameters(x)
    keyval <- keyword(x)
    sub_exprs <- exprs(x)
    sub_exprs <- sub_exprs[goodCellIDs, ]
    newx <- flowFrame(exprs = sub_exprs, parameters = params, description = keyval)
  }
  return(list(anoms = fr_autoqc$anoms, frequencies = fr_frequences,
              FRnewFCS = newx,
              goodCellIDs = goodCellIDs, badCellIDs = badCellIDs,
              res_fr_QC = data.frame(second_fraction = second_fraction,
                num_obs = fr_autoqc$num_obs, period_obs = period_obs,
                period =period, badPerc = badPerc)))
}


# Plot frequency values for a list y, containing the outputs from
# the function flow_rate_check
flow_rate_plot <- function(FlowRateQC) {

  second_fraction <- FlowRateQC$res_fr_QC$second_fraction
  num_obs = FlowRateQC$res_fr_QC$num_obs
  period_obs = FlowRateQC$res_fr_QC$period_obs
  period = FlowRateQC$res_fr_QC$period
  frequencies = as.data.frame(FlowRateQC$frequencies)
  anoms = as.data.frame(FlowRateQC$anoms)
  anoms_points = as.data.frame(cbind(sec_anom = frequencies$secbin[anoms$index], count_anom = anoms$anoms))


    lines_period <- seq(1, num_obs, period_obs)
    lines_long_period <- lines_period[seq(1, length(lines_period),
                          period)]

    xgraph <- ggplot(frequencies, aes_string(x="secbin", y="tbCounts")) +
              theme_bw() + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text=element_text(size = 14)) + geom_line(colour = "red" )
   # xgraph <- xgraph + scale_x_continuous(breaks=frequencies$minbin[lines_at],
   #           expand=c(0,0))
    xgraph <- xgraph + geom_vline(xintercept=frequencies$secbin[lines_long_period], color="gray60") + geom_vline(xintercept=frequencies$secbin[lines_period], linetype = 2, color="gray60")
    xgraph <- xgraph + labs(x= "Seconds", y= paste0("Number of cells per 1/",
              1  /second_fraction, " of a second"), title= "Flow Rate")

  # Add anoms to the plot as circles.
    if(!is.null(anoms_points)){
      xgraph <- xgraph + geom_point(data=anoms_points, aes_string(x= "sec_anom", y= "count_anom"), color = "green4", size = 3, shape = 1)
    }

    return(xgraph)
}
