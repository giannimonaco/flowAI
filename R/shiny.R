#' Interactive quality control of Flow Cytometry Data
#'
#' The call of the \code{flow_iQC} function opens a Shiny application that allows 
#' to perfom a complete and interactive quality control of an FCS file. The 
#' framework of the interactive quality control is complementary to the 
#' automatic one of the \code{flow_auto_qc} function. Hence, the anomalies 
#' are manually selected from the evaluation of three main properties of 
#' flow cytometry: 1) flow rate, 2) signal acquisition, 3) dynamic range.
#' 
#' @author Chen Hao, Gianni Monaco
#' @import RColorBrewer
#' @export
#' @examples if (interactive()) flowAI::flow_iQC()
flow_iQC = function() {
    shiny::runApp(system.file('shiny', package = 'flowAI'))
}
