## remove last slash if present
strip.sep <- function(name) {
    ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,
        substr(name,1,nchar(name)-1),name)
}

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

## create new flowFrame with the parameter indicating good and bad cells
addQC <- function(QCvector, remove_from, sub_exprs, params, keyval){
    
    rs <- attr(sub_exprs, "ranges")
    rs <- c(rs, rs[1])
    sub_exprs <- cbind(sub_exprs, QCvector)
    attr(sub_exprs, "ranges") <- rs
    NN <- as.numeric(keyval["$PAR"]) + 1
    names(dimnames(sub_exprs)[[2]])[dim(sub_exprs)[2]] <- paste0("$P",NN, "N")
    pnr <- paste0("$P", NN, "R")
    pnb <- paste0("$P", NN, "B")
    pne <- paste0("$P", NN, "E")
    pnn <- paste0("$P", NN, "N")
    pns <- paste0("$P", NN, "S")
    flowCorePnRmax <- paste0("flowCore_$P", NN, "Rmax")
    flowCorePnRmin <- paste0("flowCore_$P", NN, "Rmin")
    o <- params@data
    o[length(o[,1]) + 1,] <- list(paste0("remove_from_", remove_from), "QC", 20000, 0, 20000)
    rownames(o)[length(o[,1])] <- paste("$P", NN, sep = "")
    
    outFCS <- new("flowFrame", exprs=sub_exprs, parameters=new("AnnotatedDataFrame",o), description=keyval)
    keyword(outFCS)[pnr] <- "20000"
    keyword(outFCS)[pnb] <- "32"   # keyword(outFCS)$`$P1B`
    keyword(outFCS)[pne] <- "0,0"
    keyword(outFCS)[pnn] <- paste0("remove_from_", remove_from)
    keyword(outFCS)[pns] <- "QC"
    keyword(outFCS)$`$PAR` <- NN
    keyword(outFCS)[flowCorePnRmax] <- 20000
    keyword(outFCS)[flowCorePnRmin] <- 0
    outFCS
}  

