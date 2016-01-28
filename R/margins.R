# The events on the upper margin and the outlier in the negative
# range of values are detected and removed.
#
flow_margin_check <- function(x,  margin_channels = NULL,
                      side = "both") {

  if (is.null(margin_channels)) {
    teCh <- grep("Time|time|Event|event", colnames(x), value = TRUE)
    parms <- setdiff(colnames(x), teCh)
  } else {
    if (!all(margin_channels %in% colnames(x)))
      stop("Invalid channel(s)")
    parms <- margin_channels
  }
  scatter_parms <- grep("FSC|SSC", parms, value = TRUE)

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
      min_value <- -3.5 * mad(x[neg]) / 0.6745 + median(x[neg])  # -3.5 is the default threshold
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
  badPerc <- round(length(badCellIDs)/lenx, 2)

  if (badPerc > 0.5) {
    warning(paste0(100 * badPerc, "%", " of cells are detected as margin events \n"))
  } else {
    cat(paste0(100 * badPerc, "%", " of cells are detected as margin events \n"))
  }

  params <- parameters(x)
  keyval <- keyword(x)
  sub_exprs <- exprs(x)
  sub_exprs <- sub_exprs[goodCellIDs, ]
  newx <- flowFrame(exprs = sub_exprs, parameters = params,
          description = keyval)

  return(list(FMnewFCS = newx, goodCellIDs = goodCellIDs,
              bad_lowerIDs = bad_lowerIDs, bad_upperIDs = bad_upperIDs,
              margin_events = summary_bad_cells, badPerc = badPerc,
              events = lenx))
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