# map_export
#' Convert a consensus genetic map created with mapfuser to JoinMap format
#'
#' Writes the consensus map to the JoinMap ".map" format
#' @param MF.obj The mapfuser object with filled results slot
#' @param file A connection or a string giving the file name to write the map
#' to.
#' @author Dennis van Muijen
#' @importFrom utils write.table
#' @examples
#' \dontshow{
#' tmp <- tempfile()
#' MF.obj <- list()
#' class(MF.obj) <- "mapfuser"
#' MF.obj$result$consensus_map <- data.frame(Marker = LETTERS[1:5], LG = c(1,1,1,2,2), Position = 1:5)
#' map_export(MF.obj, file = tmp)
#' file.remove(paste0(tmp,".map"))
#' }
#' \dontrun{
#' map_export(MF.obj = MF.obj, file = "consensus")
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
  mapvec
  for (i in seq_along(grp)) {
    mapvec <- rbind.data.frame(mapvec, paste(""))
    mapvec <- rbind.data.frame(mapvec, paste("group ", grp[i]), stringsAsFactors = F, sep = "")
    for (j in which(map$LG == grp[i]))
      mapvec <- rbind.data.frame(mapvec, paste(map$Marker[j], map$Position[j], " ;"), stringsAsFactors = F)
  }
  if (!is.null(file)){
    write.table(mapvec, paste0(file,".map"), sep = "", quote = FALSE, col.names = FALSE, row.names = FALSE)
  } else {
    cat("Filename not supplied..")
    return(mapvec)
  }
}
