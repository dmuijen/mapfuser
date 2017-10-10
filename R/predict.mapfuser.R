# predict.mapfuser
#' Predict centiMorgan positions from fitted gam models on the mapfuser object
#'
#' Takes a fitted thin plate regression spline produced by gam() and produces predictions of centi Morgan positions using known physical genome positions and fitted gam models.
#' @param object A mapfuser object with fitted gam models
#' @param to_predict A csv file with columns marker, chromosome ID and position in mega base pairs.
#' @param ... ingored in function call
#' @return A data frame with columns "Marker", "Chr", "Position", and "Position_physical"
#' @author Dennis van Muijen
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "Bur", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' ref_file <- list.files(fpath, pattern = "reference", full.names = TRUE)
#' MF.obj <- read_ref(MF.obj = MF.obj, ref_file = ref_file, sep = ",", header = TRUE, type = "delim")
#' MF.obj$config$chromosomes <- "1"
#' MF.obj <- genphys_fit(MF.obj, map = "Col-0_Bur-0.csv", type = "map")
#' fpath <- system.file("extdata", package="mapfuser")
#' prd <- read.table(paste0(fpath, "/BaySha_physical.csv"), sep = ",", header = TRUE)
#' MF.obj <- predict(MF.obj, prd)
#' }
#' \dontrun{
#' # Read a table with positions to interpolate and/or extrapolate
#' fpath <- system.file("extdata", package="mapfuser")
#' to_predict <- read.table(paste0(fpath, "/BaySha_physical.csv"), sep = ",", header = TRUE)
#' MF.obj <- predict(MF.obj, to_predict)
#' # Write to csv
#' write.table(MF.obj$predictions, file = "preds.csv", sep = ",", col.names = TRUE, row.names = FALSE)
#' }
#' @export

predict.mapfuser <- function(object, to_predict, ...) {
  if (class(MF.obj) != "mapfuser") {
    stop("Object should be of class mapfuser")
  }
  if (is.null(MF.obj$pspline)) {
    stop("Fit gam models first using genphys_fit()")
  }
  chromosomes <- MF.obj$config$chromosomes
  preds <- NULL
  names(to_predict) <- c("Marker", "Chr", "Position_physical")
  fitted_values <- NULL
  for (i in seq_along(chromosomes)) {
    df <- filter(to_predict, (!!quote(Chr)) == chromosomes[i])
    df$Predictions <- predict(MF.obj$pspline$gam_models[[chromosomes[i]]], df)
    fitted_values <- rbind(fitted_values, df)
  }
  fitted_values <- arrange(fitted_values, !! as.name("Chr"), !! as.name("Position_physical"))
  fitted_values$Position <- round(fitted_values$Position, digits = 2)
  MF.obj$predictions <- fitted_values
  return(MF.obj)
}
