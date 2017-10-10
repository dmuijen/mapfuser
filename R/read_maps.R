# read_maps
#' Read genetic maps
#'
#' Reads genetic maps in either delimited or JoinMap format. Maps are loaded to a mapfuser object.
#' @param mapfiles List of filenames to read
#' @param sep The field separator character, see read.delim
#' @param header A logical value indicating whether the files contains the names of variable int the first line, see read.delim()
#' @param na.strings How to interpret missing values? See read.delim()
#' @param type Type of input data, either delimited file with columns "Marker", "LG", and "Position" or a JoinMap ".map" file.
#' @param mapweights numeric vector mapweights
#' @return A mapfuser object with the genetic maps listed under the raw_data slot. Chromosome identifiers are added to the config slot
#' @author Dennis van Muijen
#' @import ggplot2 mgcv lazyeval igraph visNetwork
#' @importFrom stringi stri_split stri_split_regex
#' @importFrom plotly ggplotly
#' @importFrom stats cor na.omit predict residuals sd
#' @importFrom utils read.delim sessionInfo tail
#' @importFrom dplyr select filter group_by distinct arrange mutate summarize n inner_join bind_rows slice
#' @examples
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "Col", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE,
#' mapweights = rep(1,7), type = "delim")
#' @export

read_maps <- function(mapfiles = NULL, sep = NULL, header = TRUE, na.strings = "NA", type = c("delim", "JoinMap"),
                      mapweights = NULL) {
  if (type == "delim" & is.null(sep)) {
    stop("specify file separator for type = 'delim'")
  }
  '.' <- NULL
  MF.obj <- vector("list", 6)
  class(MF.obj) <- "mapfuser"
  names(MF.obj) <- c("raw_data", "QC", "result", "config", "ref_map", "pspline")

  if (type == "delim") {
    for (i in seq_along(mapfiles)) {
      MF.obj$raw_data[[i]] <- read.delim(mapfiles[i], sep = sep, header = header, stringsAsFactors = FALSE, na.strings = na.strings)
      names(MF.obj$raw_data[[i]]) <- c("Marker", "LG", "Position")
      if (is.null(mapweights)) {
        MF.obj$raw_data[[i]]$Weight <- 1
      } else {
        MF.obj$raw_data[[i]]$Weight <- mapweights[i]
      }
    }
  }
  if (type == "JoinMap") {
    for (j in seq_along(mapfiles)) {
      MF.obj$raw_data[[j]] <- read_joinmap(mapfiles[i])
      names(MF.obj$raw_data[[j]]) <- c("Marker", "LG", "Position")
      if (is.null(mapweights)) {
        MF.obj$raw_data[[j]]$Weight <- 1
      } else {
        MF.obj$raw_data[[j]]$Weight <- mapweights[i]
      }
    }
  }
  if (!is.null(mapweights)) {
    if ((mapweights %>%
      length()) != (mapfiles %>%
      length())) {
      stop("Number of maps and mapweights do not correspond")
    }
  }
  for (k in seq_along(mapfiles)) {
    if ((MF.obj$raw_data[[k]]$Position %>%
      check_cM() %>%
      length()) != 0) {
      stop(paste0("Non-numeric centiMorgan position(s) specified in map ", mapfiles[k]))
    }
  }
  MF.obj$config$chromosomes <- (MF.obj$raw_data %>%
    bind_rows() %>%
    select("LG"))$LG %>%
    gsub("\\..*", "", .) %>%
    unique()
  mapnames <- lapply(mapfiles, function(x) stri_split(x, fixed = "/") %>%
    unlist() %>%
    tail(n = 1)) %>%
    unlist()
  names(MF.obj$raw_data) <- mapnames
  MF.obj$config$mapnames <- mapnames
  MF.obj$raw_data <- lapply(MF.obj$raw_data, function(x) {
    x$Marker <- as.character(x$Marker)
    x$LG <- as.character(x$LG)
    x$Position <- as.numeric(x$Position) 
    return(x)
  })
  MF.obj$input <- MF.obj$raw_data
  names(MF.obj$input) <- mapnames
  MF.obj$config$session <- sessionInfo()
  return(MF.obj)
}
