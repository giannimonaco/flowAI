#' Automatic quality control on FCS files.
#'
#' For a set of FCS files, \emph{flow_auto_qc} performs the quality control in three steps.
#' It automatically detects and removes anomalies in the flow rate, in the
#' signal acquisition and in the dynamic range.
#'
#' @param fcsfiles It can be a character vector with the filenames of the FCS files,
#' a flowSet or a flowFrame.
#' @param remove_from Decide which steps have to be considered for the removal of
#' anomalies. The full quality control is performed by default with the option \code{"all"}.
#' Other options are: \code{'FR_FS' | 'FR_FM' | 'FS_FM' | 'FR' | 'FS' | 'FM'}.
#' The abbreviation \emph{FR} commands the QC for flow rate,
#' \emph{FS} commands the QC for signal acquisition, and \emph{FM} commands
#' the QC for the margins of the dynamic range.
#' @param timeCh Character string corresponding to the name of the Time Channel
#' in the set of FCS files. By default is \code{NULL} and the name is retrieved
#' automatically.
#' @param second_fractionFR The fraction of a second that is used to split
#' the time channel in order to recreate the flow rate. Set it to 
#' \code{"timestep"} if you wish to recreate the flow rate at the maximum 
#' resolution allowed by the flow cytometry instruments. Usuallythe timestep 
#' corresponds to 0.01, however, to shorten the running time of the analysis the 
#' fraction used by default is 0.1, corresponding to 1/10 of a second.
#' @param alphaFR The level of statistical significance used to
#' accept anomalies detected by the ESD method. The default value is \code{0.01}.
#' @param decompFR Logical indicating whether the flow rate should be decomposed 
#' in the trend and cyclic components. Default is \code{TRUE} and the ESD 
#' outlier detection will be done on the trend component penalized with the 
#' magnitude of the cyclic component. If it is \code{FALSE} the ESD outlier 
#' detection will be done on the original flow rate.
#' @param ChRemoveFS Add a character vector with the names or
#' the pattern of the names of the channels that you do not want to include
#' in the signal quality control. The default option excludes the 
#' scatter parameters using \code{c("FSC", "SSC")}. If you want
#' to include all the  parameters in the analysis use \code{NULL}.
#' @param outlierFS logical indicating whether to remove outliers from 
#' the intensity values
#' before performing the changepoint analysis. The default is \code{FALSE}.
#' @param pen_valueFS The value of the penalty for the changepoint detection
#' algorithm. This can be a numeric value or text giving the formula to use,
#' for instance it is possible to use the character string "1.5*log(n)",
#' where n indicates the number of cells in the FCS file. The higher the
#' penalty value the less strict is the detection of the anomalies.
#' The default value is \code{200}.
#' @param max_cptFS The maximum number of changepoints that can be detected
#' for each channel. The default value is \code{3}.
#' @param ChFM A character vector that indicates which channels
#' need to undergo quality control of the dynamic range. The default
#' option is \code{NULL} and with it all the channels are selected for the analysis.
#' @param sideFM Select if the checking of the dynamic range should be performed in the
#' lower side, upper side or both sides of the acquitision range. Choose between:
#' \code{'both' | 'upper' | 'lower'}. The default value is \code{both}.
#' @param neg_valuesFM Scalar indicating the method to use for the removal of the 
#' anomalies from the lower limit of the dynamic range. Use \code{1} to remove 
#' negative outliers or use \code{2} to truncate the negative values to the cut-off
#' indicated in the FCS file. 
#' @param html_report Character string that will be added to the filename of the 
#' fcs analyzed to name a html document reporting the results of the quality control. 
#' The default is \code{"_QC"}. If you do not want to generate a report 
#' use \code{FALSE}.
#' @param mini_report create a text file with the percentage of anomalies detected in 
#' each FCS file analysed. The default is \code{"_QCmini"}. If you prefer not to generate 
#' the mini report use \code{FALSE}.
#' @param fcs_highQ Character string that will be added to the filename of the 
#' fcs analyzed to name the new FCS file that will be generated containing only the 
#' cells that passed the quality control. The default is \code{"_HighQ"}. If you do not 
#' want to generate a new FCS file use \code{FALSE}.
#' @param fcs_lowQ Character string that will be added to the filename of the 
#' fcs analyzed to name the new FCS file that will be generated containing only the 
#' cells that do not passed the quality control. By default the FCS file is not 
#' generated and the argument is set as \code{FALSE} 
#' @param folder_results Character string for the name of the directory that will contain 
#' the results. The default is \code{"resultsQC"}. If you intend return the results 
#' in the main directory use \code{FALSE}.
#' @return It will return the FCS files containing only the cells that
#' passed the quality checks. By default it will return a report indicating the cells
#' that were removed in the flow rate, signal acquisition over time and anomalous marginal
#' events on the dynamic range. It will also return a new FCS file containing
#' only the cells that did not pass the quality check.
#' @author Gianni Monaco 
#' @examples
#'
#' ## a sample file
#' data(Bcells)
#'
#' ## quality control
#' flow_auto_qc(Bcells[[1]], html_report = FALSE, mini_report = FALSE, fcs_highQ = FALSE, folder_results = FALSE)
#'
#' @import flowCore
#' @import ggplot2
#' @import plyr
#' @importFrom changepoint cpt.meanvar
#' @importFrom scales pretty_breaks
#' @import knitr
#' @import reshape2
#' @export
flow_auto_qc <- function(fcsfiles, remove_from = "all",
     timeCh = NULL, second_fractionFR = 0.1, alphaFR = 0.01, decompFR = TRUE, 
     ChRemoveFS = c("FSC", "SSC"), outlierFS = FALSE, pen_valueFS = 200, 
     max_cptFS = 3, ChFM = NULL, sideFM = "both", neg_valuesFM = 1, 
     html_report = "_QC", mini_report = "QCmini", fcs_highQ = "_HighQ", 
     fcs_lowQ = FALSE, folder_results = "resultsQC") {
  ## load the data
  if( is.character(fcsfiles) ){
    set <- read.flowSet(files = fcsfiles)
  }else if( class(fcsfiles) == "flowSet"){
    set <- fcsfiles
  }else if( class(fcsfiles) == "flowFrame" ){
    set <- as(fcsfiles,"flowSet")
  }else{
   stop("Use as first argument a flowSet or a character vector with the name of the fcs files to analyse")
  }

  N_cell_set <- flow_set_qc(set)
  area.color <- rep("red", length(set))

  FSbinSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
  if (missing(timeCh) || is.null(timeCh)) {
    timeCh <- findTimeChannel(set[[1]])
  }
  if (length(nchar(timeCh)) == 0 || is.null(timeCh)) {
    warning("Impossible to retreive it automatically. The quality control can be performed only on signal and dynamic range.")
  }

  # in some cases, especially if the FCS file has been modified, there
  # could be more than one slots for the Timestep parameter.  the first in
  # numerical order should be the original value.
  word <- which(grepl("TIMESTEP", names(set[[1]]@description),
                ignore.case = TRUE))
  timestep <- as.numeric(set[[1]]@description[[word[1]]])
  if( length(timestep) == 0 ){
    warning("Timestep object did not found in the FCS file and it was set to 0.01. Graphs labels indicating time might not be correct", call. =FALSE)
    timestep <- 0.01
  }
  if( second_fractionFR == "timestep" ){
      second_fractionFR <- timestep
  }else if( second_fractionFR < timestep ){
      stop("The argument second_fractionFR must be greater or equal to timestep.", call. =FALSE)
  }
  
  if(folder_results != FALSE){
      dir.create(folder_results, showWarnings = FALSE)
     # setwd(file.path(getwd(), folder_results))
  }
  
  for (i in 1:length(set)) {

    filename_ext <- description(set[[i]])$FILENAME
    filename <- sub("^([^.]*).*", "\\1", filename_ext)
    if (html_report != FALSE) {
        reportfile <- paste0(getwd(), .Platform$file.sep, 
            ifelse(folder_results != FALSE, paste0(folder_results, .Platform$file.sep), ""),
            filename, html_report, ".html")
    }
    if (mini_report != FALSE) {
        minireport <- paste0(getwd(), .Platform$file.sep, 
            ifelse(folder_results != FALSE, paste0(folder_results, .Platform$file.sep), ""),
            mini_report, ".txt")
        if(!file.exists(minireport)){
            write.table(t(c("Name file", "n. of events", "% anomalies", "analysis from",
                "% anomalies flow Rate",  "% anomalies Signal",  "% anomalies Margins")),
                minireport, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
        }
    }
    if (fcs_highQ != FALSE) {
        good.fcs.file <- paste0(getwd(), .Platform$file.sep,
            ifelse(folder_results != FALSE, paste0(folder_results, .Platform$file.sep), ""),
            filename, fcs_highQ, ".fcs")
    }
    if (fcs_lowQ != FALSE) {
        bad.fcs.file <- paste0(getwd(), .Platform$file.sep,
            ifelse(folder_results != FALSE, paste0(folder_results, .Platform$file.sep), ""),
            filename, fcs_lowQ, ".fcs")
    }

    # select different color for the analyzed FCS in the set plot
    area <- area.color
    area[i] <- "blue"

    # check if the FCS file is ordered by time
    if (length(nchar(timeCh)) != 0 && !is.null(timeCh)) {
        ordFCS <- ord_fcs_time(set[[i]], timeCh)
    }else{ ordFCS <- set[[i]] }
   # ordFCS <- set[[i]]   ## TEMPORARY
    origin_cellIDs <- 1:nrow(ordFCS)

    ### Describe here the arguments for the functions of the flow Rate and Flow Signal
    FR_bin_arg <- list( second_fraction = second_fractionFR, timeCh = timeCh, 
                  timestep = timestep)
    FR_QC_arg <- list( alpha = alphaFR, use_decomp = decompFR)
    FS_bin_arg <- list( binSize = FSbinSize, timeCh = timeCh, timestep = timestep )
    FS_QC_arg <- list( ChannelRemove = ChRemoveFS, pen_valueFS, max_cptFS, outlierFS )
    FM_QC_arg <- list( margin_channels = ChFM , side= sideFM, neg_values = neg_valuesFM)

    #### The actual analysis is performed here
    if (length(nchar(timeCh)) != 0 && !is.null(timeCh)) {
      FlowRateData <- do.call(flow_rate_bin, c(ordFCS, FR_bin_arg ))
      FlowRateQC <- do.call(flow_rate_check, c(ordFCS, list(FlowRateData), FR_QC_arg ))
      } else if(remove_from == "all"){
          remove_from <-"FS_FM" 
      } else if(remove_from == "FR_FS"){
          remove_from <-"FS"
      } else if(remove_from == "FR_FM"){
          remove_from <-"FM"
      } else if( remove_from == "FR"){
          stop("it is not possible to perfrom the flow rate check without a time channel")
      }
      
      FlowSignalData <- do.call(flow_signal_bin, c(ordFCS,FS_bin_arg))
      FlowSignalQC <- do.call(flow_signal_check, c(ordFCS,list(FlowSignalData),FS_QC_arg))
      FlowMarginQC <- do.call(flow_margin_check, c(ordFCS, FM_QC_arg))

      ####selection of good cells
      if(remove_from == "all"){
       goodCellIDs <- intersect(FlowRateQC$goodCellIDs, intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs))
       analysis <- "Flow Rate, Flow Signal and Flow Margin"
      }else if(remove_from == "FR_FS"){
        goodCellIDs <- intersect(FlowRateQC$goodCellIDs, FlowSignalQC$goodCellIDs)
        analysis <- "Flow Rate and Flow Signal"
      }else if(remove_from == "FR_FM"){
        goodCellIDs <- intersect(FlowRateQC$goodCellIDs, FlowMarginQC$goodCellIDs)
        analysis <- "Flow Rate and Flow Margin"
      }else if(remove_from == "FS_FM"){
        goodCellIDs <- intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs)
        analysis <- "Flow Signal and Flow Margin"
      }else if(remove_from == "FR"){
        goodCellIDs <- FlowRateQC$goodCellIDs
        analysis <- "Flow Rate"
      }else if(remove_from == "FS"){
        goodCellIDs <- FlowSignalQC$goodCellIDs
        analysis <- "Flow Signal"
      }else if(remove_from == "FM"){
        goodCellIDs <- FlowMarginQC$goodCellIDs
        analysis <- "Flow Margin"
      }

      badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)
      totalBadPerc <- round(length(badCellIDs)/length(origin_cellIDs), 2)
      if (fcs_highQ != FALSE || fcs_lowQ != FALSE) {
        params <- parameters(ordFCS)
        keyval <- keyword(ordFCS)
        sub_exprs <- exprs(ordFCS)
      }
      if (length(badCellIDs) > 1 & fcs_highQ != FALSE) {
      good_sub_exprs <- sub_exprs[goodCellIDs, ]
      goodfcs <- flowFrame(exprs = good_sub_exprs,
        parameters = params, description = keyval)
    suppressWarnings(write.FCS(goodfcs, good.fcs.file))
    }
    if (length(badCellIDs) > 1 & fcs_lowQ != FALSE) {
      bad_sub_exprs <- sub_exprs[badCellIDs, ]
      badfcs <- flowFrame(exprs = bad_sub_exprs,
        parameters = params,description = keyval)
      suppressWarnings(write.FCS(badfcs, bad.fcs.file))
    }
     # write.table(as.vector(badCellIDs), file= paste0("resultsQC/", filename, "_bad.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)  # TEMPORARY
    if (length(nchar(timeCh)) == 0 || is.null(timeCh)) {
        FlowRateQC<-list()
        FlowRateQC$res_fr_QC$badPerc <- 0
    }
    if (mini_report != FALSE) {
        write.table(t(c(filename, as.integer(dim(set[[i]])[1]), 
        totalBadPerc * 100, analysis, FlowRateQC$res_fr_QC$badPerc * 100, 
        FlowSignalQC$Perc_bad_cells$badPerc_tot * 100,  
        FlowMarginQC$badPerc * 100)), minireport, sep="\t", 
        append=TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
     }
     if (html_report != FALSE) {
       h_FS_graph <- round(0.4 * (ncol(ordFCS) - length(ChRemoveFS)),1)
       if (!is.null(ChRemoveFS)){ 
           ChannelRemovedFS <- as.character(grep(paste(ChRemoveFS, collapse="|"),
               ordFCS@parameters$name, value = TRUE))
       }
       template_path <- system.file("rmd","autoQC_report.Rmd", package='flowAI')
       knit2html(template_path, output = reportfile, force_v1 = TRUE)
     }
  }
}
