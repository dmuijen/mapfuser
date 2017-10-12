# map_qc
#' Wrapper function of genetic map cleaning per linkage group.
#'
#' The raw data is first splitted to separate linkage groups if sublinkage groups exist (e.g LG 1.1 and 1.2).
#' Subsequently a graph is created from the adjacency matrix that counts the number of overlapping markers between the set of genetic maps.
#' Calculations are performed for each chromosome separately. Taken quality control steps are printed to the console and can be visualised using the plot function.
#' @param MF.obj A mapfuser object genetics maps loaded and optionally a reference map
#' @param anchors Number of minimum overlapping anchors marker between at least one other genetic map. At least 3 are required.
#' @return The input object is returned with filled QC slot containing genetic maps after quality control. Used parameters and inverted or names of removed data are saved to the config slot.
#' @author Dennis van Muijen, Ram Kumar Basnet
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "-1", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- map_qc(MF.obj)
#' }
#' \dontrun{
#' MF.obj <- map_qc(MF.obj = MF.obj, anchors = 3)
#' #Graphical overview of how different genetic maps are connected by overlapping markers
#' plot(MF.obj, which = "mapnetwork", chr = 1) ## Multiple chromosomes not supported
#' ## A minimal spanning tree using the number of anchors as edge weight,
#' plot(MF.obj, which = "mst", chr = 1)
#' #Visualize inverted maps
#' plot(MF.obj, which = "genetic_maps", maps = c("Col-0_Cvi-0.csv","Col-0_Sha.csv"), chr = 1:3)
#' }
#' @export

map_qc <- function(MF.obj, anchors = 3) {
  MF.obj$config$n.anchors <- anchors
  MF.obj$config$cohesion <- FALSE
  if (class(MF.obj) != "mapfuser") {
    stop("Object should be of class mapfuser")
  }
  if (is.null(MF.obj$raw_data)) {
    stop("Load map(s) to mapfuser object first using read.maps")
  }
  MF.obj$config$LPmerge_ready <- FALSE
  # Split maps to (sub) linkage group for LPmerge integration
  MF.obj <- map_split(MF.obj)
  # Check for sufficient anchors to orient and integrate maps
  MF.obj <- check_anchors(MF.obj, anchors = anchors)

  # Check if there are sets of maps that cannot be integrated
  MF.obj <- map_cohesion(MF.obj)

  # Correct map orientations
  if (MF.obj$config$cohesion) {
    MF.obj <- map_orient(MF.obj)

    MF.obj$config$LPmerge_ready <- TRUE
    cat("QC Finished")
  }
  return(MF.obj)
}
