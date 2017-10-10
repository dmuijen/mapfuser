# map_cohesion
#' Graph check whether maps can be integrated
#' @param MF.obj A mapfuser object
#' @return The input mapfuser object if check OK, error otherwise
#' @author Dennis van Muijen

map_cohesion <- function(MF.obj) {
  MF.obj$QC$maps <- MF.obj$QC$maps[lapply(MF.obj$QC$maps, function(x) dim(x)[1]) > 0]
  if (is.null(MF.obj$QC$qcnetworks)) {
    stop("No networks to calculate cohesion")
  }
  MF.obj$QC$mst <- lapply(MF.obj$QC$qcnetworks, minimum.spanning.tree)
  map_sets <- lapply(MF.obj$QC$mst, cohesion) %>%
    unlist()
  if (any(map_sets < 1)) {
    bad_lg <- names(map_sets)[map_sets < 1] %>%
      as.numeric() %>%
      sort() %>%
      as.character()
    cat("Non integratable sets of maps found for chromosome(s).. ", paste(bad_lg, collapse = ", "), "\n")
    cat("Plot minimum spanning tree for details", "\n")
    cat("Remove a set of unconnected linkage groups to continue or start analysis from scratch excluding non-integratable sets", "\n")
    return(MF.obj)
  }
  MF.obj$config$cohesion <- TRUE
  return(MF.obj)
}
