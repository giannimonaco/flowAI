# Guess which channel captures time in a exprs, flowFrame or flowset
findTimeChannel <- function(xx) {
    time <- grep("^Time$", colnames(xx), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(time)) {
        if (is(xx, "flowSet") || is(xx, "ncdfFlowList"))
            xx <- exprs(xx[[1]]) else if (is(xx, "flowFrame"))
                xx <- exprs(xx)
            cont <- apply(xx, 2, function(y) all(sign(diff(y)) >= 0))
            time <- names(which(cont))
    }
    if (!length(time) || length(time) > 1)
        time <- NULL
    return(time)
}

# Check if the Fcs file is ordered according to time otherwise it order it.
ord_fcs_time <- function(x, timeCh){

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


