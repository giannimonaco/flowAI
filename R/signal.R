# A flowFrame object is splitted in bins with equal number of events
# and for each bin the median is calculated.
#
# @param x: a flowFrame object
# @param channels: channel names for the selected markers
# @param binSize: the size of bin
flow_signal_bin <- function(x, channels = NULL, binSize = 500,
  timeCh = timeCh, timestep = timestep, TimeChCheck = TimeChCheck) {

  ## some sanity checking
  if (!is(x, "flowFrame"))
    stop("'x' needs to be of class 'flowFrame'")

  if (is.null(channels) || missing(channels) || is.na(channels)) {
    parms <- setdiff(colnames(x), timeCh)
  } else {
    if (!all(channels %in% colnames(x)))
      stop("Invalid channel(s)")
    parms <- channels
  }

  ### Retriving time and expression info
  exp <- exprs(x)
  if (!is.null(TimeChCheck)) {
      timex <- seq(from = 0, length.out = nrow(x), by = 0.1)
  }else{
    timex <- exp[, timeCh]
  }
  yy <- exp[, parms]  # channels data
  idx <- c(1:nrow(x))
  seconds <- timex * timestep
  lenSec <- length(seconds)  # num of time ticks
  uniSeconds <- unique(seconds)  # num of unique time tick
  nrBins <- floor(lenSec/binSize)  # num of bins

  cf <- c(rep(1:nrBins, each = binSize), rep(nrBins + 1, lenSec - nrBins * binSize))  # id bins
  stopifnot(length(cf) == lenSec)
  tmpx <- split(seconds, cf)
  xx2 <- sapply(tmpx, mean)  # mean of each time bin  (x axis)
  yy2 <- as.matrix(ddply(as.data.frame(yy), .(cf), colwise(median)))[, -1]

  return(list(exprsBin = cbind(timeSec = xx2, yy2), cellBinID = data.frame(cellID = idx, binID = cf),
    bins = length(unique(cf)), binSize = binSize))
}

#########################################################################
# Detection of shifts in the median intensity signal detected
# by the laser of the flow cytometry over time
flow_signal_check <- function(x, FlowSignalData, ChannelExclude = NULL,
  pen_valueFS = pen_valueFS, maxSegmentFS = maxSegmentFS, outlier_remove = FALSE) {

  fs_cellBinID <- FlowSignalData$cellBinID
  fs_res <- FlowSignalData$exprsBin
  ### log transformation.
  # fs_res[which(fs_res <= 1 & fs_res >= -1)] <- 0
  # fs_res[which(fs_res > 1)] <- log(fs_res[which(fs_res > 1)])
  # fs_res[which(fs_res < -1)] <- -log(abs(fs_res[which(fs_res < -1)]))

  teCh <- grep("Time|time|TIME|Event|event|EVENT", colnames(fs_res), value = TRUE)
  parms <- setdiff(colnames(fs_res), teCh)

  ##### scale and sum the value of each channel and find the outliers
  scale_sign <- apply(fs_res[, parms],2, scale)
  sum_sign <- apply(abs(scale_sign),1, sum)
  outup_tr <- (+3.5 * mad(sum_sign) + (0.6745 * median(sum_sign)))/0.6745
  FS_out <- which(sum_sign >outup_tr)

  #### Remove channel from the changepoint analysis
  if (!is.null(ChannelExclude)) {
    ChannelExclude_COMP <- grep(paste(ChannelExclude, collapse="|"),
                      colnames(fs_res), value = TRUE)
  #  cat(paste0("The channels whose signal acquisition will not be checked are: ",
   #   paste(ChannelExclude_COMP, collapse = ", "), ". \n"))
    parms <- setdiff(parms, ChannelExclude_COMP)
  }

  if(outlier_remove){
    ##transform outliers to median values
    if( length(FS_out) != 0 ){
      fs_res_adj <- apply(fs_res[, parms], 2, function(x) {
        med <- median(x)
        x[FS_out] <- med
        return(x)
      })
     # badPerc_out <- round((length(FS_out)/nrow(fs_res)),4)
     # cat(paste0(badPerc_out * 100, "% of outliers found in channels' signal. \n"))
    }else{
      fs_res_adj <- fs_res[,parms]
    #  cat("0% of outliers found in channels' signal. \n")
     # badPerc_out <- 0
    }

  cpt_res <- suppressWarnings(cpt.meanvar(t(fs_res_adj),
             pen.value = pen_valueFS, Q = maxSegmentFS,
             penalty =  "Manual" , test.stat = "Normal",
             method = "BinSeg", param.estimates = FALSE))
  }else{

  cpt_res <- suppressWarnings(cpt.meanvar(t(fs_res[, parms]),
      pen.value = pen_valueFS, Q = maxSegmentFS,
      penalty =  "Manual" , test.stat = "Normal",
      method = "BinSeg", param.estimates = FALSE))
  badPerc_out <- 0
  }

  list_seg <- lapply(1:length(cpt_res), function(x) {
    # it retrieves the changepoints that has been detected for each channel
    cpt_res[[x]]@cpts
  })
  names(list_seg) <- parms
  list_seg <- as.list(list_seg)

  list_cpt <- union(1, sort(unique(unlist(list_seg[parms]))))  # retrieve and order all the changepoints
  diff_cpt <- sapply(2:length(list_cpt), function(n) {
    # calculate the difference between each segment
    list_cpt[n] - list_cpt[n - 1]
  })
  max_diff <- which(diff_cpt == max(diff_cpt))
  max_seg <- c(list_cpt[max_diff], list_cpt[max_diff + 1] )  # selecting the biggest segment

  list_seg <- lapply(list_seg, function(x) setdiff(x, nrow(fs_res)))

  len_cpt <- sapply(list_seg, length)
  nam_cpt <- gsub("<|>","",names(len_cpt))
  nozero_cpt <- as.numeric(which(len_cpt != 0))
  zero_cpt <- as.numeric(which(len_cpt == 0))
  if(length(nozero_cpt) == 0){
    ch_no_cpt <- nam_cpt[zero_cpt]
    tab_cpt <- NULL
  }else{
   # cat(paste("Changepoint(s) detected in the channels: ",
   # paste(names(len_cpt[nozero_cpt]), collapse = ", "), sep = ""), fill = TRUE)
    ch_cpt <- list_seg[nozero_cpt]
    ch_no_cpt <- nam_cpt[zero_cpt]

    max_n_cpt <- max(sapply(ch_cpt, length))
    tab_cpt <- ldply(ch_cpt, function(x) c(x, rep(NA, max_n_cpt - length(x))),
      .id = NULL)
    rownames(tab_cpt) <- nam_cpt[nozero_cpt]
    tab_cpt <- as.matrix(tab_cpt)
    tab_cpt[which(is.na(tab_cpt))] <- ""
    colnames(tab_cpt) <- 1:length(tab_cpt[1, ])
  }
  # percentage bad cell detected with the changepoint method
  # badPerc_cp <- round(1 - ((max_seg[2] - max_seg[1])/(length(fs_res[, 1]) - 1)),4)
  
  # retrieve ID of good cells
  if(outlier_remove){
    fs_cellBinID2 <- fs_cellBinID[which(!fs_cellBinID[, 2] %in% FS_out),]
    badPerc_out <- round(1 - nrow(fs_cellBinID2)/nrow(fs_cellBinID),4)
    goodCellIDs <- fs_cellBinID2[which(fs_cellBinID2[, 2] >= max_seg[1] &
                                         fs_cellBinID2[, 2] <= max_seg[2]), 1]
    badPerc_tot <- round(1 - length(goodCellIDs)/nrow(fs_cellBinID),4)
  }else{
    goodCellIDs <- fs_cellBinID[which(fs_cellBinID[, 2] >= max_seg[1] &
                                        fs_cellBinID[, 2] <= max_seg[2]), 1]
    badPerc_tot <- round(1 - length(goodCellIDs)/nrow(fs_cellBinID),4)
  }
  
  cat(paste0(100 * badPerc_tot, "% of anomalous cells detected in signal acquisition check. \n"))

  params <- parameters(x)
  keyval <- keyword(x)
  sub_exprs <- exprs(x)
  check_goodcellsID(goodCellIDs)
  sub_exprs <- sub_exprs[goodCellIDs, , drop = FALSE]  ## check if the Id Correspond!
  newx <- flowFrame(exprs = sub_exprs, parameters = params, description = keyval)

  return(list(FSnewFCS = newx, exprsBin = FlowSignalData$exprsBin, Perc_bad_cells = data.frame(badPerc_tot,badPerc_out),
        goodCellIDs = goodCellIDs, tab_cpt = tab_cpt, ch_no_cpt =ch_no_cpt,
        segm = max_seg, FS_out = FS_out, outlier_remove = outlier_remove))
}


# Plot the flourescence intensity for each channel of a flowFrame
# over time, highlighting the wider segment that do not show shifts
# of the median intensity
# @param exprsBin give the exprsBin object from FlowSignalData
# @param segm give the segm object from FlowSignalQC
# @param FS_out give the FS_out object from FlowSignalQC
flow_signal_plot <- function(FlowSignalQC) {

  exprsBin <- FlowSignalQC$exprsBin
  segm <- FlowSignalQC$segm
  FS_out <- FlowSignalQC$FS_out
  outlier_remove <- FlowSignalQC$outlier_remove

  binID <- 1:nrow(exprsBin)
  teCh <- grep("Time|time|Event|event", colnames(exprsBin), value = TRUE)
  parms <- setdiff(colnames(exprsBin), teCh)
  dataORIG <- exprsBin[, parms]  # first channel is time

  if(length(FS_out) != 0){
    data <- apply(dataORIG, 2, function(x){
     overMAX <- FS_out[x[FS_out] > max(x[-FS_out])]
     overMIN <- FS_out[x[FS_out] < min(x[-FS_out])]
     x[overMAX] <- max(x[-FS_out])
     x[overMIN] <- min(x[-FS_out])
     return(x)
    })
    data <- as.data.frame(data)
  }else{
    data <- as.data.frame(dataORIG)
  }
  data$binID <- binID
  longdata <- melt(data, id.vars = "binID", variable.name = "marker",
    value.name = "value")
  FS_graph <- ggplot(longdata, aes_string(x = "binID", y = "value", col = "marker"),
      environment = environment()) + labs(x = "Bin ID",
      y = "Median Intensity value") +
      facet_grid(marker ~ ., scales = "free") + theme_bw() +
      theme(strip.text.y = element_text(angle = 0,
      hjust = 1), axis.text.y = element_text(size = 6),
      legend.position = "none") +
      scale_x_continuous(breaks= pretty_breaks(n =10)) +
      scale_y_continuous(breaks= pretty_breaks(n =3)) +
      geom_rect(aes(xmin = segm[1], xmax = segm[2], ymin = -Inf,
      ymax = Inf), fill = "khaki1", linetype = 0) + geom_line()
  # Add anoms to the plot as circles.
  if(outlier_remove){
    longdata_out <- melt(data[FS_out,], id.vars = "binID",
      variable.name = "marker",value.name = "value")
    FS_graph <- FS_graph + geom_point(data=longdata_out, aes_string(x= "binID",
                y= "value"), color = "green4", size = 2, shape = 1)
  }
  return(FS_graph)
}
