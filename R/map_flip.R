# map_flip
#' Invert a linkage group
#'
#' @param x Linkage group to invert
#' @return  Linkage group with inverted genetic positions
#' @author Dennis van Muijen
#' @export
map_flip <- function(x) {
  the_max <- max(x$Position)
  x$Position <- (x$Position - the_max) * -1
  x <- arrange(x, !! as.name("Position"))
  return(x)
}
