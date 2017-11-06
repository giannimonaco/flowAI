## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ---- message=FALSE------------------------------------------------------
require(flowAI)

## ---- collapse = TRUE----------------------------------------------------
data(Bcells)
Bcells

## ---- eval=FALSE---------------------------------------------------------
#  setwd(...)
#  fcsfiles <- dir(".", pattern="*fcs$")

## ---- eval=FALSE---------------------------------------------------------
#  flow_auto_qc(Bcells)  # using a flowSet
#  flow_auto_qc(Bcells[[1]]) # using a flowFrame
#  flow_auto_qc(fcsfiles) # using a character vector

## ---- eval=FALSE---------------------------------------------------------
#  GbLimit <- 2    # decide the limit in gigabyte for your batches of FCS files
#  size_fcs <- file.size(fcsfiles)/1024/1024/1024    # it calculates the size in gigabytes for each FCS file
#  groups <- ceiling(sum(size_fcs)/GbLimit)
#  cums <- cumsum(size_fcs)
#  batches <- cut(cums, groups)

## ---- eval=FALSE---------------------------------------------------------
#  for(i in 1:groups){
#      flow_auto_qc(fcsfiles[which(batches == levels(batches)[i])])
#  }

