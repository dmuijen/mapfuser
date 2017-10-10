# remove_manual
#' Manually remove a linkage group within a specific map and update man networks and minimum spanning tree
#'
#' The quality control passed list of genetic maps may be manually curated further with remove_manual() when the map merging process identified a linkage group witihn a map that gives a high Root Mean Square Error.
#' In the case of interspecific crosses one complete map could be better to exclude all together.
#' @param MF.obj A mapfuser object genetics maps after map_orient() has been performed
#' @param to_remove A list of linkage group ID's to remove
#' @return The input object is returned in which the linkage group ID have been removed. Manually removed linkage group ID are saved to the removed_LGIDs config slot
#' @author Dennis van Muijen
#' @examples
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "-1", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- map_qc(MF.obj)
#' ## Remove two linkage groups manually
#' to_remove <- c("Col-0_Blh-1.csv_1","Col-0_Blh-1.csv_2" )
#' MF.obj <- remove_manual(MF.obj, to_remove)
#' @export

remove_manual <- function(MF.obj, to_remove) {
  for (i in seq_along(to_remove)) {
    MF.obj$QC$maps[[to_remove[i]]] <- NULL
  }
  MF.obj$input <- split(
    MF.obj$QC$maps %>%
      bind_rows() %>%
      select("Marker", "LG", "Position", "Weight"),
    (MF.obj$QC$maps %>%
      bind_rows())$id
  )
  MF.obj$config$manually_removed <- c(MF.obj$config$manually_removed, to_remove)
  return(MF.obj)
}
