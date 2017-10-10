# map_export
#' Convert a consensus genetic map created with mapfuser to JoinMap format
#'
#' Writes the consensus map to the JoinMap ".map" format
#' @param MF.obj The mapfuser object with filled results slot
#' @param file Path to the output file
#' @author Dennis van Muijen
#' @importFrom utils write.table
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "-1", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- map_qc(MF.obj)
#' MF.obj$config$chr <- "1"
#' MF.obj <- LPmerge_par(MF.obj = MF.obj, n.cores = 1,
#' max.interval = 1, max.int_sel = "auto", weights = NULL)
#' file <- paste(tempdir(), "/consensus.map", sep="")
#' map_export(MF.obj, file)
#' }
#' \dontrun{
#' ## Read maps
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "Col", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE,
#' mapweights = rep(1,7), type = "delim")
#' 
#' ## Run map_qc
#' MF.obj <- map_qc(MF.obj, anchors = 3)
#' 
#' ## Construct consensus map
#' MF.obj <- LPmerge_par(MF.obj = MF.obj, n.cores = 2,
#' max.interval = 1:3, max.int_sel = "auto", weights = NULL)
#' 
#' ## Export to JoinMap format
#' file <- paste(tempdir(), "/consensus.map", sep="")
#' map_export(MF.obj, file)
#' }
#' @export

map_export <- function(MF.obj, file = NULL) {
  if (class(MF.obj) != "mapfuser") {
    stop("MF.obj should be of class mapfuser")
  }
  map <- MF.obj$result$consensus_map
  grp <- unique(map$LG)
  mapvec <- paste("; ", date(), collapse = "")
  mapvec <- rbind.data.frame(mapvec, paste(
    "; ngrp = ", length(unique(map$LG)),
    ", nloc = ", nrow(map), sep = ""
  ), stringsAsFactors = F)
  for (i in seq_along(grp)) {
    mapvec <- rbind.data.frame(mapvec, paste(""))
    mapvec <- rbind.data.frame(mapvec, paste("group ", grp[i]), stringsAsFactors = F, sep = "")
    for (j in which(map$LG == grp[i]))
      mapvec <- rbind.data.frame(mapvec, paste(map$Marker[j], map$Position[j], " ;"), stringsAsFactors = F)
  }
  if (!is.null(file)){
    write.table(mapvec, file, sep = "", quote = FALSE, col.names = FALSE, row.names = FALSE)
  } else {
    cat("Filename not supplied.. returning genetic map in JoinMap format")
    return(mapvec)
  }
}
