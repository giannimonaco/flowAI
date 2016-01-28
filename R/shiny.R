#' A Shiny app to perform interactive quality control on Flow Cytometry Data
#'
#' This function allows to interactively perform the quality control of .fcs files in a Shiny
#' app. The quality control framework of the Shiny app is equivalent to the
#' one of the \code{flow_auto_qc} function. It is composed of three main steps that allow
#' to check anomalies in the 1) flow rate, 2) signal acquisition, 3) dynamic range.
#' @import RColorBrewer
#' @export
#' @examples if (interactive()) flowAI::flow_iQC()
flow_iQC = function() {
    shiny::runApp(system.file('shiny', package = 'flowAI'))
}
