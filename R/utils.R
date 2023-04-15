check_goodcellsID <- function(x) {
  if (length(x) == 0) {
    warning("There are not good events", call. = FALSE)
  }
}


check_badCellIDs <- function(x) {
  if (length(x) == 0) {
    warning("There are not bad events", call. = FALSE)
  }
}
