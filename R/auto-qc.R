#' Automatic quality control on .fcs files.
#'
#' For a set of .fcs files, \emph{flow_auto_qc} performs the quality control in three steps.
#' It automatically detects and removes anomalies in the flow rate, in the
#' signal acquisition and in the dynamic range.
#'
#' @param fcsfiles It can be a character vector with the filenames of the .fcs files,
#' a flowSet or a flowFrame.
#' @param remove_from Decide which steps have to be considered for the removal of
#' anomalies. The full quality control is performed by default with the option \code{"all"}.
#' Other options are: \code{'FR_FS' | 'FR_FM' | 'FS_FM' | 'FR' | 'FS' | 'FM'}.
#' The abbreviation \emph{FR} commands the QC for flow rate,
#' \emph{FS} commands the QC for signal acquisition, and \emph{FM} commands
#' the QC for the margins of the dynamic range.
#' @param highQ_suffix Character string to add to the original filenames
#' to name the new .fcs files containing the cells that passed the quality control.
#' The default is \code{NULL}, that automatically adds the suffix \emph{_highQ}
#' to the original name.
#' @param report logical indicating whether the report has to be released.
#' The default is \code{TRUE}.
#' @param report_suffix Character string to adds to the original filenames
#' to name the reports that will be released by the quality control.
#' The default is \code{NULL}, that automatically add the suffix \emph{_reportQC}
#' to the original name.
#' @param save_highQ logical indicating whether the high quality .fcs files have to be
#' saved. The default is \code{TRUE}.
#' @param save_lowQ logical indicating whether the low quality .fcs files have to be
#' saved. The default is \code{TRUE}.
#' @param lowQ_suffix Character string to add to the original filenames
#' to name the new .fcs files containing the cells that did not pass the quality control.
#' The default is \code{NULL}, that automatically adds the suffix \emph{_LowQ}
#' to the original name.
#' @param timeCh Character string corresponding to the name of the Time Channel
#' in the set of .fcs files. By default is \code{NULL} and the name is retrieved
#' automatically.
#' @param second_fractionFR The fraction of a second that is used to split
#' the time channel to recreate the flow rate. The fraction used by default is 1/10 of a second.
#' @param long_periodFR An integer indicating the number of sections the flow rate
#' should be split into and indipendently analysed. The default setting is \code{1}
#' and the flow rate is analyzed in its entirety.
#' @param periodFR An integer that indicates the number of segments a long period
#' in the flow rate has to be split into in order for the anomaly to be detected.
#' The default value is \code{30}.
#' @param alphaFR The level of statistical significance used to
#' accept anomalies. The default value is \code{0.05}.
#' @param ChRemoveFS Add a character vector with the names or
#' the pattern of the names of the channels that you do not want to include
#' in the signal quality control. The default option is \code{NULL}. If you want
#' to exclude the scatter parameters from the analysis use \code{c("FCS", "SSC")}.
#' @param outlierFS logical indicating whether to remove outliers from the intensity values
#' before performing the changepoint analysis. The default is \code{FALSE}.
#' @param pen_valueFS The value of the penalty for the changepoint detection
#' algorithm. This can be a numeric value or text giving the formula to use,
#' for instance it is possible to use the character string "1.5*log(n)",
#' where n indicates the number of cells in the .fcs file. The higher the
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
#' @return It will return the .fcs files containing only the cells that
#' passed the quality checks. By default it will return a report indicating the cells
#' that were removed in the flow rate, signal acquisition over time and anomalous marginal
#' events on the dynamic range. It will also return a new .fcs file containing
#' only the cells that did not pass the quality check.
#' @author Gianni Monaco with substantial contribution from Chen Hao
#' @examples
#'
#' ## a sample file
#' data(Bcells)
#'
#' ## quality control
#' flow_auto_qc(Bcells[[1]], report = FALSE,  save_highQ = FALSE, save_lowQ = FALSE)
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
     highQ_suffix = NULL, report = TRUE, report_suffix = NULL,
     save_highQ = TRUE, save_lowQ = TRUE, lowQ_suffix = NULL,
     timeCh = NULL, second_fractionFR = 0.1, long_periodFR = 1,
     periodFR = 30, alphaFR = 0.05, ChRemoveFS = NULL, outlierFS = FALSE,
     pen_valueFS = 200, max_cptFS = 3,
     ChFM = NULL, sideFM = "both") {
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
    stop("Time channel needed. Impossible to retreive it automatically. Check fcs file.")
  }

  # in some cases, especially if the FCS file has been modified, there
  # could be more the one slots for the Timestep parameter.  the first in
  # numerical order should be the original value.
  word <- which(grepl("TIMESTEP", names(set[[1]]@description),
                ignore.case = TRUE))
  timestep <- as.numeric(set[[1]]@description[[word[1]]])
  if( length(timestep) == 0 ){
    warning("Timestep object did not found in the FCS file and it was setted to 0.1. graphs labels might not correspond to reality.", call. =FALSE)
    timestep <- 0.01
  }


  for (i in 1:length(set)) {

    filename_ext <- description(set[[i]])$FILENAME
    filename <- sub("^([^.]*).*", "\\1", filename_ext)
    if (is.null(highQ_suffix)) {
      good.fcs.file <- paste0(getwd(), .Platform$file.sep,
                               filename, "_HighQ.fcs")
    }else{
      good.fcs.file <- paste0(getwd(), .Platform$file.sep,
        filename, highQ_suffix, ".fcs")
    }
    if (is.null(lowQ_suffix)) {
      bad.fcs.file <- paste0(getwd(), .Platform$file.sep,
                         filename, "_LowQ.fcs")
    }else{
      good.fcs.file <- paste0(getwd(), .Platform$file.sep,
        filename, lowQ_suffix, ".fcs")
    }
    if (is.null(report_suffix)){
      reportfile <- paste0(getwd(), .Platform$file.sep,
                                filename, "_reportQC.html")
    }else{
      reportfile <- paste0(getwd(), .Platform$file.sep,
        filename, report_suffix, ".html")
    }

    # select different color for the analyzed FCS in the set plot
    area <- area.color
    area[i] <- "blue"

    # check if the FCS file is ordered by time
    ordFCS <- ord_fcs_time(set[[i]], timeCh)
    origin_cellIDs <- 1:nrow(ordFCS)

    ### Describe here the arguments for the functions of the flow Rate and Flow Signal
    FR_bin_arg <- list( second_fraction = second_fractionFR,
      timeCh = timeCh, timestep = timestep)
    FR_QC_arg <- list( long_period = long_periodFR, period = periodFR,
      alpha = alphaFR )
    FS_bin_arg <- list( binSize = FSbinSize, timeCh = timeCh,
                  timestep = timestep )
    FS_QC_arg <- list( ChannelRemove = ChRemoveFS, pen_valueFS, max_cptFS, outlierFS )
    FM_QC_arg <- list( margin_channels = ChFM , side= sideFM)

    #### The actual analysis is performed here
      FlowRateData <- do.call(flow_rate_bin, c(ordFCS, FR_bin_arg ))
      FlowRateQC <- do.call(flow_rate_check, c(ordFCS, list(FlowRateData), FR_QC_arg ))
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
      if (length(badCellIDs) > 1 & save_highQ == TRUE) {
      params <- parameters(ordFCS)
      keyval <- keyword(ordFCS)
      sub_exprs <- exprs(ordFCS)

      good_sub_exprs <- sub_exprs[goodCellIDs, ]
      goodfcs <- flowFrame(exprs = good_sub_exprs,
        parameters = params, description = keyval)
    suppressWarnings(write.FCS(goodfcs, good.fcs.file))
    }
    if (length(badCellIDs) > 1 & save_lowQ == TRUE) {
      bad_sub_exprs <- sub_exprs[badCellIDs, ]
      badfcs <- flowFrame(exprs = bad_sub_exprs,
        parameters = params,description = keyval)
      suppressWarnings(write.FCS(badfcs, bad.fcs.file))
    }

    if (report == TRUE) {
       h_FS_graph <- round(0.4 * (ncol(ordFCS) - length(ChRemoveFS)),1)
       template_path <- system.file('rmd/autoQC_report.Rmd',
       package='flowAutoQC')

       knit2html(template_path, output = reportfile)
    }
  }
}
