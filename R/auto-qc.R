#' Automatic quality control of flow cytometry data.
#'
#' For a set of FCS files, \emph{flow_auto_qc} performs a complete and automatic 
#' quality control. It consists in the detection and removal of
#' anomalies by checking three properties of flow cytometry: 1) flow rate, 
#' 2) signal acquisition, 3) dynamic range.
#'
#' @param fcsfiles It can be a character vector with the filenames of the FCS files,
#' a flowSet or a flowFrame.
#' @param remove_from Select from which of the three steps the anomalies have to be 
#' excluded in the high quality FCS file. The default option \code{"all"} removes the 
#' anomalies from all the three steps. Alternatively, you can use: 
#' \code{"FR_FS", "FR_FM", "FS_FM", "FR", "FS", "FM"}, to remove the anomalies only 
#' on a subset of the steps where \emph{FR} stands for the flow rate, \emph{FS} stands
#' for signal acquisition and \emph{FM} stands for dynamic range.
#' @param timeCh Character string corresponding to the name of the Time Channel
#' in the set of FCS files. By default is \code{NULL} and the name is retrieved
#' automatically.
#' @param second_fractionFR The fraction of a second that is used to split
#' the time channel in order to recreate the flow rate. Set it to 
#' \code{"timestep"} if you wish to recreate the flow rate at the maximum 
#' resolution allowed by the flow cytometry instrument. Usually, the timestep 
#' corresponds to 0.01, however, to shorten the running time of the analysis the 
#' fraction used by default is 0.1, corresponding to 1/10 of a second.
#' @param alphaFR The level of statistical significance used to
#' accept anomalies detected by the ESD method. The default value is \code{0.01}.
#' @param decompFR Logical indicating whether the flow rate should be decomposed 
#' in the trend and cyclical components. Default is \code{TRUE} and the ESD 
#' outlier detection will be executed on the trend component penalized by the 
#' magnitude of the cyclical component. If it is \code{FALSE} the ESD outlier 
#' detection will be executed on the original flow rate.
#' @param ChRemoveFS Add a character vector with the names or name portions 
#' of the channels that you want to exclude from the signal acquisition check. 
#' The default option, \code{c("FSC", "SSC")}, excludes the scatter parameters. 
#' If you want to include all the parameters in the analysis use \code{NULL}.
#' @param outlierFS logical indicating whether outliers have to be removed
#' before the changepoint detection of the signal acquisition check. 
#' The default is \code{FALSE}.
#' @param pen_valueFS The value of the penalty for the changepoint detection
#' algorithm. This can be a numeric value or text giving the formula to use;
#' for instance, you can use the character string \code{"1.5*log(n)"},
#' where n indicates the number of cells in the FCS file. The higher the
#' penalty value the less strict is the detection of the anomalies.
#' The default is \code{200}.
#' @param max_cptFS The maximum number of changepoints that can be detected
#' for each channel. The default is \code{3}.
#' @param ChFM A character vector that indicates which channels need to include
#' for the dynamic range check. The default option is \code{NULL} and 
#' with it all the channels are selected for the analysis.
#' @param sideFM Select whether the dynamic range check has to be executed on  
#' both limits, the upper limit or the lower limit. Use one of the options:
#' \code{"both", "upper", "lower"}. The default is \code{"both"}.
#' @param neg_valuesFM Scalar indicating the method to use for the removal of the 
#' anomalies from the lower limit of the dynamic range. Use \code{1} to remove 
#' negative outliers or use \code{2} to truncate the negative values to the cut-off
#' indicated in the FCS file. 
#' @param html_report Suffix to be added to the FCS filename to name the HTML 
#' report of the quality control. The default is \code{"_QC"}. If you do not 
#' want to generate a report use \code{FALSE}.
#' @param mini_report Suffix to be added for the filename of the TXT 
#' file containing the percentage of anomalies detected in each FCS file 
#' analyzed. The default is \code{"_QCmini"}. If you do not want to generate 
#' the mini report use \code{FALSE}.
#' @param fcs_QC Suffix to be added for the filename of the new FCS 
#' containing a new channel where the low quality events have a random value 
#' between 10,000 to 20,000 (as for flowClean). The default is 
#' \code{"_QC"}. If you do not want to generate the high quality FCS file 
#' use \code{FALSE}.
#' @param fcs_highQ Suffix to be added for the filename of the new FCS 
#' containing only the events that passed the quality control. The default 
#' is \code{FALSE} and hence the high quality FCS file is not generated.
#' @param fcs_lowQ Suffix to be added for the filename of the new FCS 
#' containing only the events that did not pass the quality control. The default 
#' is \code{FALSE} and hence the low quality FCS file is not generated.
#' @param folder_results Character string used to name the directory that 
#' contains the results. The default is \code{"resultsQC"}. If you intend 
#' to return the results in the main directory use \code{FALSE}.
#' @return A complete quality control is performed on flow cytometry data in FCS
#' format. By default the analysis returns a directory named \code{resultsQC}
#' containing:
#' 1. a set of new FCS files with a new parameter to gate out the low quality events
#' 2. a set of HTML reports, one for each FCS file, that include graphs and table 
#' indicating where the anomalies were detected,
#' 3. a single TXT file reporting the percentage of events removed in each FCS file.
#' 
#' @author Gianni Monaco, Chen Hao 
#' @examples
#'
#' ## a sample dataset as flowSet object
#' data(Bcells)
#'
#' ## quality control on a flowFrame object
#' flow_auto_qc(Bcells[[1]], html_report = FALSE, mini_report = FALSE, fcs_QC = FALSE, folder_results = FALSE)
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
     html_report = "_QC", mini_report = "QCmini", fcs_QC = "_QC", fcs_highQ = FALSE, 
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
  if (is.null(timeCh)) {
    warning("Impossible to retrieve the time channel automatically. The quality control can only be performed on signal acquisition and dynamic range.", call. =FALSE)
  }

  # in some cases, especially if the FCS file has been modified, there
  # could be more than one slots for the Timestep parameter. the first in
  # numerical order should be the original value.
  word <- which(grepl("TIMESTEP", names(set[[1]]@description),
                ignore.case = TRUE))
  timestep <- as.numeric(set[[1]]@description[[word[1]]])
  if( !length(timestep) ){
    warning("The TIMESTEP keyword was not found and hence it was set to 0.01. Graphs labels indicating time might not be correct", call. =FALSE)
    timestep <- 0.01
  }
  if( second_fractionFR == "timestep" ){
      second_fractionFR <- timestep
  }else if( second_fractionFR < timestep ){
      stop("The argument second_fractionFR must be greater or equal to timestep.", call. =FALSE)
  }
  
  if(folder_results != FALSE){
      dir.create(folder_results, showWarnings = FALSE)
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
    if (fcs_QC != FALSE) {
        QC.fcs.file <- paste0(getwd(), .Platform$file.sep,
            ifelse(folder_results != FALSE, paste0(folder_results, .Platform$file.sep), ""),
            filename, fcs_QC, ".fcs")
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

    # check the time channel of the file
    if (!is.null(timeCh)) {
        if (length(unique(exprs(set[[i]])[, timeCh])) == 1){
            cat("The time channel contain a single value. It cannot be used to recreate the flow rate. \n")
            warning(paste0("The time channel in ", filename_ext, " contains a single value. It cannot be used to recreate the flow rate. \n"), call. =FALSE)
            TimeChCheck <- "single_value"
        }else{
            TimeChCheck <- NULL 
        }
    }else{
        TimeChCheck <- "NoTime"
    }
    
    # order FCS file is a proper Time channel is present
    if (is.null(TimeChCheck)) {
      ordFCS <- ord_fcs_time(set[[i]], timeCh)
    }else{ 
      ordFCS <- set[[i]] 
    }
    
    origin_cellIDs <- 1:nrow(ordFCS)

    ### Describe here the arguments for the functions of the flow Rate and Flow Signal
    FR_bin_arg <- list( second_fraction = second_fractionFR, timeCh = timeCh, 
                  timestep = timestep)
    FR_QC_arg <- list( alpha = alphaFR, use_decomp = decompFR)
    FS_bin_arg <- list( binSize = FSbinSize, timeCh = timeCh, timestep = timestep, TimeChCheck = TimeChCheck)
    FS_QC_arg <- list( ChannelRemove = ChRemoveFS, pen_valueFS, max_cptFS, outlierFS )
    FM_QC_arg <- list( margin_channels = ChFM , side= sideFM, neg_values = neg_valuesFM)

    #### The actual analysis is performed here
      if (is.null(TimeChCheck)) {
        FlowRateData <- do.call(flow_rate_bin, c(ordFCS, FR_bin_arg ))
        FlowRateQC <- do.call(flow_rate_check, c(ordFCS, list(FlowRateData), FR_QC_arg ))
      }else{
        FlowRateQC<-list()
        FlowRateQC$goodCellIDs <- origin_cellIDs
        FlowRateQC$res_fr_QC$badPerc <- 0
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
      totalBadPerc <- round(length(badCellIDs)/length(origin_cellIDs), 4)
      if (fcs_QC != FALSE || fcs_highQ != FALSE || fcs_lowQ != FALSE) {
        sub_exprs <- exprs(ordFCS)
        params <- parameters(ordFCS)
        keyval <- keyword(ordFCS)
      }
      if (fcs_QC != FALSE ){
          QCvector <- FlowSignalData$cellBinID[,"binID"]
          QCvector[badCellIDs] <- runif(length(badCellIDs), min=10000, max=20000) 
          newFCS <- addQC(QCvector, remove_from, sub_exprs, params, keyval)
          suppressWarnings(write.FCS(newFCS, QC.fcs.file))  
      }
      if (length(badCellIDs) > 0 & fcs_highQ != FALSE) {
      good_sub_exprs <- sub_exprs[goodCellIDs, ]
      goodfcs <- flowFrame(exprs = good_sub_exprs,
        parameters = params, description = keyval)
    suppressWarnings(write.FCS(goodfcs, good.fcs.file))
    }
    if (length(badCellIDs) > 0 & fcs_lowQ != FALSE) {
      bad_sub_exprs <- sub_exprs[badCellIDs, ]
      badfcs <- flowFrame(exprs = bad_sub_exprs,
        parameters = params,description = keyval)
      suppressWarnings(write.FCS(badfcs, bad.fcs.file))
    }
    if (mini_report != FALSE) {
        write.table(t(c(filename, as.integer(dim(set[[i]])[1]), 
        totalBadPerc * 100, analysis, FlowRateQC$res_fr_QC$badPerc * 100, 
        FlowSignalQC$Perc_bad_cells$badPerc_tot * 100,  
        FlowMarginQC$badPerc * 100)), minireport, sep="\t", 
        append=TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
     }
     if (html_report != FALSE) {
       h_FS_graph <- round(0.4 * (ncol(ordFCS)),1)
       if (!is.null(ChRemoveFS)){ 
           ChannelRemovedFS <- as.character(grep(paste(ChRemoveFS, collapse="|"),
               ordFCS@parameters$name, value = TRUE))
       }
       template_path <- system.file("rmd","autoQC_report.Rmd", package='flowAI')
       knit2html(template_path, output = reportfile, force_v1 = TRUE)
     }
  }
}
