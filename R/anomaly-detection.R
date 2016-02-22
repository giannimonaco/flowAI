# Detects anomalies in a time series using Cyclic hybrid ESD (C-H-ESD).
#
# Anomaly Detection Using Cyclic Hybrid ESD Test ----GM
#
# A technique for detecting anomalies in univariate time
# series where the input is a series of observations.
# @name AnomalyDetection
# @param x Time series as a column data frame, list, or vector,
#  where the column consists of the observations.
# @param max_anoms Maximum number of anomalies that C-H-ESD will
# detect as a percentage of the data. The value can be from 0 to 1.
# @param direction Directionality of the anomalies to be detected.
# Options are: \code{'pos' | 'neg' | 'both'}.
# @param alpha The level of statistical significance with which
# to accept or reject anomalies.
# @param use_decomp If set to \code{'FALSE'} it gives the possibility 
# to detect outliers with the generalized ESD method on the orginal data.  
# By default is set to \code{'TRUE'} and time series decomposition
# is performed before the analysis.
# @param period Defines the number of observations in a single
# period, and used during seasonal decomposition.
# @param e_value Add an additional column to the anoms output
# containing the expected value.
# @param verbose Additionally printing for debugging.

# @details
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, decomposition components, and
# optionally expected values.}
# @return One can save \code{anoms} to a file in the following fashion:
# \code{write.csv(<return list name>[["anoms"]], file=<filename>)}
# @references Rosner, B., (May 1983), "Percentage Points for a
# Generalized ESD Many-Outlier Procedure", Technometrics, 25(2),
# pp. 165-172.
# @export


anomaly_detection = function(x, max_anoms=0.49, direction='both', alpha=0.01, use_decomp = TRUE, period=1, verbose = FALSE){
    
    # Check for supported inputs types
    if(is.vector(x) && is.numeric(x)) {
        x <- ts(x, frequency = period)
    } else if(is.ts(x)) {
    } else {
        stop("data must be a time series object or a vector that holds numeric values.")
    }
    
    # Handle NAs
    if (length(rle(is.na(c(NA,x,NA)))$values)>3){
        stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
    } else {
        x <- na.omit(x)
    }
    
    # Sanity check all input parameterss
    if(max_anoms > .49){
        stop(paste("max_anoms must be less than 50% of the data points (max_anoms =", round(max_anoms*length(x), 0), " data_points =", length(x),")."))
    }
    if(!direction %in% c('pos', 'neg', 'both')){
        stop("direction options are: pos | neg | both.")
    }
    if(!(0.01 <= alpha || alpha <= 0.1)){
        print("Warning: alpha is the statistical significance level, and is usually between 0.01 and 0.1")
    }
    if(is.null(period)){
        stop("Period must be set to the number of data points in a single period")
    }

    ############## -- Main analysis: Perform C-H-ESD -- #################
    # -- Step 1: Decompose data. This will return two more components: trend and cycle   
    if(use_decomp){
        x_cf <- cffilter(x)
        #med_t <- trunc(median(x_cf$trend))
        med_t <- trunc(median(x))
        sign_n <- sign(x_cf$trend - med_t)
        sign_n[which(sign_n == 0)] <-1
        # add the absolute values of the cycle component to the absolute values of the centered trend component. The signs are then added again 
        x_2 <- as.vector(trunc(abs(x - med_t) + abs(x_cf$cycle)) * sign_n)
    } else {
        x_2 <- as.vector(x - median(x))
    }
    
    anomaly_direction = switch(direction,
        "pos" = data.frame(one_tail=TRUE, upper_tail=TRUE), # upper-tail only (positive going anomalies)
        "neg" = data.frame(one_tail=TRUE, upper_tail=FALSE), # lower-tail only (negative going anomalies)
        "both" = data.frame(one_tail=FALSE, upper_tail=TRUE)) # Both tails. Tail direction is not actually used.
    

    n <- length(x_2)
    data_det <- data.frame(index = 1:n, values = x_2, or_values = x)
    # Maximum number of outliers that C-H-ESD can detect (e.g. 49% of data)
    max_outliers <- trunc(n*max_anoms)
    func_ma <- match.fun(median)
    func_sigma <- match.fun(mad)
    R_idx <- 1L:max_outliers
    num_anoms <- 0L
    one_tail <- anomaly_direction$one_tail
    upper_tail <- anomaly_direction$upper_tail
    # Compute test statistic until r=max_outliers values have been
    # removed from the sample.
    for (i in 1L:max_outliers){
        if(verbose) message(paste(i,"/", max_outliers,"completed"))
        
        if(one_tail){
            if(upper_tail){
                ares <- data_det[[2L]] - func_ma(data_det[[2L]])
            } else {
                ares <- func_ma(data_det[[2L]]) - data_det[[2L]]
            }
        } else {
            ares = abs(data_det[[2L]] - func_ma(data_det[[2L]]))
        }
        
        # protect against constant time series
        data_sigma <- func_sigma(data_det[[3L]])   
        # the standard deviation has to be calculated from the orginal 
        # distribution because otherwise it would be affected to match 
        # by the cycle component
        if(data_sigma == 0) 
            break
        
        ares <- ares/data_sigma
        R <- max(ares)
        
        temp_max_idx <- which(ares == R)[1L]
        
        R_idx[i] <- data_det[[1L]][temp_max_idx]
        
        data_det <- data_det[-which(data_det[[1L]] == R_idx[i]), ]
        
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
        all_data <- data.frame(index = 1:n, anoms = x)
        anoms_data <- subset(all_data, (all_data[[1]] %in% R_idx))
    } else {
        anoms_data <- NULL
    }
    return (list(anoms = anoms_data, num_obs = n))
}

