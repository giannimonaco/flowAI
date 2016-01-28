# Detects anomalies in a time series using S-H-ESD.
#
# Args:
#	 data: Time series to perform anomaly detection on.
#	 k: Maximum number of anomalies that S-H-ESD will detect as a percentage of the data.
#	 alpha: The level of statistical significance with which to accept or reject anomalies.
#	 num_obs_per_period: Defines the number of observations in a single period, and used during seasonal decomposition.
#	 use_decomp: Use seasonal decomposition during anomaly detection.
#	 use_esd: Uses regular ESD instead of hybrid-ESD. Note hybrid-ESD is more statistically robust.
#	 one_tail: If TRUE only positive or negative going anomalies are detected depending on if upper_tail is TRUE or FALSE.
#	 upper_tail: If TRUE and one_tail is also TRUE, detect only positive going (right-tailed) anomalies. If FALSE and one_tail is TRUE, only detect negative (left-tailed) anomalies.
#	 verbose: Additionally printing for debugging.
# Returns:
#   A list containing the anomalies (anoms) and decomposition components (stl).
detect_anoms <- function(data, k = 0.49, alpha = 0.05, num_obs_per_period = NULL,
  use_decomp = TRUE, use_esd = FALSE, one_tail = TRUE,
  upper_tail = TRUE, verbose = FALSE) {


  if(is.null(num_obs_per_period)){
    stop("must supply period length for time series decomposition")
  }

  num_obs <- nrow(data)

  # Check to make sure we have at least two periods worth of data for anomaly context
  if(num_obs < num_obs_per_period * 2){
    stop("Anom detection needs at least 2 periods worth of data")
  }

  # Handle NAs
  if (length(rle(is.na(c(NA,data[[2L]],NA)))$values)>3){
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  } else {
    data <- na.omit(data)
  }

  # -- Step 1: Decompose data. This returns a univarite remainder which will be used for anomaly detection. Optionally, we might NOT decompose.
  data_decomp <- stl(ts(data[[2L]], frequency = num_obs_per_period),
    s.window = "periodic", robust = TRUE)

  # Remove the seasonal component, and the median of the data to create the univariate remainder
  data <- data.frame(timestamp = data[[1L]], count = (data[[2L]]-data_decomp$time.series[,"seasonal"]-median(data[[2L]])))

  # Store the smoothed seasonal component, plus the trend component for use in determining the "expected values" option
  data_decomp <- data.frame(timestamp=data[[1L]], count=(as.numeric(trunc(data_decomp$time.series[,"trend"]+data_decomp$time.series[,"seasonal"]))))

  # Maximum number of outliers that S-H-ESD can detect (e.g. 49% of data)
  max_outliers <- trunc(num_obs*k)

  if(max_outliers == 0){
    stop(paste0("With longterm=TRUE, AnomalyDetection splits the data into 2 week periods by default. You have ", num_obs, " observations in a period, which is too few. Set a higher piecewise_median_period_weeks."))
  }

  func_ma <- match.fun(median)
  func_sigma <- match.fun(mad)

  ## Define values and vectors.
  n <- length(data[[2L]])
  R_idx <- 1L:max_outliers

  num_anoms <- 0L

  # Compute test statistic until r=max_outliers values have been
  # removed from the sample.
  for (i in 1L:max_outliers){
    if(verbose) message(paste(i,"/", max_outliers,"completed"))

    if(one_tail){
      if(upper_tail){
        ares <- data[[2L]] - func_ma(data[[2L]])
      } else {
        ares <- func_ma(data[[2L]]) - data[[2L]]
      }
    } else {
      ares = abs(data[[2L]] - func_ma(data[[2L]]))
    }

    # protect against constant time series
    data_sigma <- func_sigma(data[[2L]])
    if(data_sigma == 0)
      break

    ares <- ares/data_sigma
    R <- max(ares)

    temp_max_idx <- which(ares == R)[1L]

    R_idx[i] <- data[[1L]][temp_max_idx]

    data <- data[-which(data[[1L]] == R_idx[i]), ]

    ## Compute critical value.
    if(one_tail){
      p <- 1 - alpha/(n-i+1)
    } else {
      p <- 1 - alpha/(2*(n-i+1))
    }

    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))

    if(R > lam)
      num_anoms <- i
  }

  if(num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
  } else {
    R_idx = NULL
  }

  return(list(anoms = R_idx, stl = data_decomp))
}


# Anomaly Detection Using Seasonal Hybrid ESD Test ----GM
#
# A technique for detecting anomalies in seasonal univariate time
# series where the input is a series of observations.
# @name AnomalyDetectionVec
# @param x Time series as a column data frame, list, or vector,
#  where the column consists of the observations.
# @param max_anoms Maximum number of anomalies that S-H-ESD will
# detect as a percentage of the data.
# @param direction Directionality of the anomalies to be detected.
# Options are:
# \code{'pos' | 'neg' | 'both'}.
# @param alpha The level of statistical significance with which
# to accept or reject anomalies.
# @param period Defines the number of observations in a single
# period, and used during seasonal decomposition.
# @param e_value Add an additional column to the anoms output
# containing the expected value.
# @param longterm_period Defines the number of observations for
# which the trend can be considered flat. The value should be an
# integer multiple of the number of observations in a single period.
# This increases anom detection efficacy for time series that are
# greater than a month.
# @details
# \code{longterm_period} This option should be set when the input
# time series is longer than a month. The option enables the
# approach described in Vallis, Hochenbaum, and Kejariwal (2014).\cr\cr
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, values, and
# optionally expected values.}
# @return One can save \code{anoms} to a file in the following fashion:
# \code{write.csv(<return list name>[["anoms"]], file=<filename>)}
# @references Vallis, O., Hochenbaum, J. and Kejariwal, A., (2014)
# "A Novel Technique for Long-Term Anomaly Detection in the Cloud",
# 6th USENIX, Philadelphia, PA.
# @references Rosner, B., (May 1983), "Percentage Points for a
# Generalized ESD Many-Outlier Procedure", Technometrics, 25(2),
# pp. 165-172.
# @export
anomaly_detection = function(x, max_anoms=0.10, direction='pos', alpha=0.05, period=NULL, e_value=TRUE, longterm_period=NULL){

    # add timestamps
   x <- data.frame(timestamp=c(1:length(x)), count=x)

    # Sanity check all input parameterss
    if(max_anoms > .49){
        stop(paste("max_anoms must be less than 50% of the data points (max_anoms =", round(max_anoms*length(x[[2]]), 0), " data_points =", length(x[[2]]),")."))
    }
    if(!direction %in% c('pos', 'neg', 'both')){
        stop("direction options are: pos | neg | both.")
    }
    if(!(0.01 <= alpha || alpha <= 0.1)){
        print("Warning: alpha is the statistical significance, and is usually between 0.01 and 0.1")
    }
    if(is.null(period)){
        stop("Period must be set to the number of data points in a single period")
    }
    if(!is.logical(e_value)){
        stop("e_value must be either TRUE (T) or FALSE (F)")
    }

    # -- Main analysis: Perform S-H-ESD

    num_obs <- length(x[[2]])

    if(max_anoms < 1/num_obs){
        max_anoms <- 1/num_obs
    }

    # -- Setup for longterm time series

    # If longterm is enabled, break the data into subset data frames and store in all_data,
    if(!is.null(longterm_period)){
        all_data <- vector(mode="list", length=ceiling(length(x[[1]])/(longterm_period)))
        # Subset x into two week chunks
        for(j in seq(1,length(x[[1]]), by=longterm_period)){
            start_index <- x[[1]][j]
            end_index <- min((start_index + longterm_period - 1), num_obs)
            # if there is at least longterm_period left, subset it, otherwise subset last_index - longterm_period
            if((end_index - start_index + 1) == longterm_period){
                all_data[[ceiling(j/(longterm_period))]] <- subset(x, x[[1]] >= start_index & x[[1]] <= end_index)
            }else{
                all_data[[ceiling(j/(longterm_period))]] <- subset(x, x[[1]] > (num_obs-longterm_period) & x[[1]] <= num_obs)
            }
        }
    }else{
        # If longterm is not enabled, then just overwrite all_data list with x as the only item
        all_data <- list(x)
    }

    # Create empty data frames to store all anoms and seasonal+trend component from decomposition
    all_anoms <- data.frame(timestamp=numeric(0), count=numeric(0))
    seasonal_plus_trend <- data.frame(timestamp=numeric(0), count=numeric(0))

    # Detect anomalies on all data (either entire data in one-pass, or in 2 week blocks if longterm=TRUE)
    for(i in 1:length(all_data)) {

        anomaly_direction = switch(direction,
                                   "pos" = data.frame(one_tail=TRUE, upper_tail=TRUE), # upper-tail only (positive going anomalies)
                                   "neg" = data.frame(one_tail=TRUE, upper_tail=FALSE), # lower-tail only (negative going anomalies)
                                   "both" = data.frame(one_tail=FALSE, upper_tail=TRUE)) # Both tails. Tail direction is not actually used.

        # detect_anoms actually performs the anomaly detection and returns the results in a list containing the anomalies
        # as well as the decomposed components of the time series for further analysis.
        s_h_esd_timestamps <- detect_anoms(all_data[[i]], k=max_anoms, alpha=alpha, num_obs_per_period=period, use_decomp=TRUE, use_esd=FALSE,
                                           one_tail=anomaly_direction$one_tail, upper_tail=anomaly_direction$upper_tail, verbose=FALSE)

        # store decomposed components in local variable and overwrite s_h_esd_timestamps to contain only the anom timestamps
        data_decomp <- s_h_esd_timestamps$stl
        s_h_esd_timestamps <- s_h_esd_timestamps$anoms

        # -- Step 3: Use detected anomaly timestamps to extract the actual anomalies (timestamp and value) from the data
        if(!is.null(s_h_esd_timestamps)){
            anoms <- subset(all_data[[i]], (all_data[[i]][[1]] %in% s_h_esd_timestamps))
        } else {
            anoms <- data.frame(timestamp=numeric(0), count=numeric(0))
        }

        all_anoms <- rbind(all_anoms, anoms)
        seasonal_plus_trend <- rbind(seasonal_plus_trend, data_decomp)
    }

    # Cleanup potential duplicates
    all_anoms <- all_anoms[!duplicated(all_anoms[[1]]), ]
    seasonal_plus_trend <- seasonal_plus_trend[!duplicated(seasonal_plus_trend[[1]]), ]

    # Calculate number of anomalies as a percentage
    anom_pct <- (length(all_anoms[[2]]) / num_obs) * 100

    # Store expected values if set by user
    if(e_value) {
      anoms <- data.frame(index=all_anoms[[1]], anoms=all_anoms[[2]], expected_value=subset(seasonal_plus_trend[[2]], seasonal_plus_trend[[1]] %in% all_anoms[[1]]))
    } else {
      anoms <- data.frame(index = all_anoms[[1]], anoms = all_anoms[[2]])
    }

       # If there are no anoms, then let's exit
    if(anom_pct == 0){
        print("No anomalies detected.")
        return (list(anoms=NULL, num_obs = num_obs, period = period))
    }
    return (list(anoms = anoms, num_obs = num_obs, period = period))
}
