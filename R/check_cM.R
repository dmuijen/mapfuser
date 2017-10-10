# check_cM
#' internal position check
#'
#' @param x vector of positions
#' @return index of non-numeric items
#' @author Dennis van Muijen

check_cM <- function(x) {
  non_numeric <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(non_numeric & !is.na(x))
}
