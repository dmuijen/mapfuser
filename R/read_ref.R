# read_ref
#' Load a reference map
#'
#' Read a reference map, either a reference genetic map or a map with physical genome positions
#' @param MF.obj A mapfuser object
#' @param sep The field separator character, see read.delim
#' @param ref_file path to reference file
#' @param header A logical value indicating whether the files contains the names of variable int the first line, see read.delim()
#' @param na.strings How to interpret missing values? See read.delim
#' @param type Type of input data, either delimited file with columns "Marker", "LG", and "Position" or a JoinMap ".map" file.
#' @return A mapfuser object with the reference map loaded to the ref_map slot. Name of the reference map is saved to the config.
#' @author Dennis van Muijen
#' @examples
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "Col", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE,
#' mapweights = rep(1,7), type = "delim")
#' ref_file <- list.files(fpath, pattern = "reference", full.names = TRUE)
#' MF.obj <- read_ref(MF.obj = MF.obj, ref_file = ref_file, sep = ",",
#' header = TRUE, na.string = NA, type = "delim")
#' @export

read_ref <- function(MF.obj = NULL, ref_file = ref_file, sep = NULL, header = TRUE, na.strings = "NA", type = c("delim", "JoinMap")) {
  if (class(MF.obj) != "mapfuser") {
    stop("MF.obj should be of class mapfuser")
  }
  if (is.null(MF.obj$raw_data)) {
    stop("Run read_maps first")
  }
  if (type == "delim") {
    MF.obj$ref_map <- read.delim(ref_file, sep = sep, header = header, stringsAsFactors = FALSE, na.strings = na.strings)
    if (ncol(MF.obj$ref_map) > 3) {
      stop("More than three columns in input - check input file")
    }
  }
  if (type == "JoinMap") {
    MF.obj$ref_map <- read_joinmap(ref_file)
  }
  names(MF.obj$ref_map) <- c("Marker", "Chr", "Position_physical")
  MF.obj$ref_map$Marker <- MF.obj$ref_map$Marker %>%
    as.character()
  MF.obj$ref_map$Chr <- MF.obj$ref_map$Chr %>%
    as.factor()
  MF.obj$config$ref_name <- stri_split_regex(ref_file, "/") %>%
    lapply(tail, n = 1) %>%
    unlist()
  return(MF.obj)
}
