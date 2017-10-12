# genphys_fit
#' Model the relationship between genetical and physical genome positions
#'
#' Fit a penalized thin plate regression spline (Wood, 2003) through genetic position in centiMorgan and physical genome positions in (Mega) base pairs.
#' @param MF.obj A mapfuser object with a reference map loaded and either a consensus map created with mapfuser or a genetic loaded with read.maps
#' @param type Fit physical genome positions vs. consensus genetic map made with mapfuser or an individual genetic map
#' @param map Name of the genetic map to use when the consensus map is not used for fitting a P-spline and recombination rate calculation
#' @param z discard invididual data points based on z-score threshold for scaled pearson residuals
#' @param chromosomes The chromosomes to fit a P-spline, default to all chromosomes
#' @return The input object is returned with added components recombination rate at a 0.1 Mbp interval,
#' the general additive model fit (gam) with penalized spline fits per chromosome,
#' and predictions of centiMorgan positions used for fitting. Markers that have been removed due to z-threshold are saved to the config slot.
#' @author Dennis van Muijen, Nathalie Dek
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "Bur", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' ref_file <- list.files(fpath, pattern = "reference", full.names = TRUE)
#' MF.obj <- read_ref(MF.obj = MF.obj, ref_file = ref_file, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- genphys_fit(MF.obj,chromosomes = 1, map = "Col-0_Bur-0.csv", type = "map")
#' }
#' \dontrun{
#' MF.obj <- genphys_fit(MF.obj, type = "consensus", z = 5, chromosomes = 1:5, map = NULL)
#' MF.obj <- genphys_fit(MF.obj, type = "map", z = 5, chromosomes = 1:5, map = "Col-0_Blh-1.csv")
#' # Plot the result
#' plot(MF.obj, which = "mareymap", maps = "consensus", chr = 1:5)
#' }
#' @export

genphys_fit <- function(MF.obj, type = c("consensus", "map"), z = 5, chromosomes = NULL, map = NULL) {
  if (class(MF.obj) != "mapfuser") {
    stop("Object should be of class mapfuser")
  }
  if (is.null(MF.obj$ref_map)) {
    stop("Load a reference map first")
  }
  if (is.null(chromosomes)) {
    chromosomes <- MF.obj$config$chromosomes
  }
  if (type == "consensus") {
    df <- inner_join(MF.obj$ref_map, MF.obj$result$consensus_map, by = "Marker") %>%
      arrange(!! as.name("Chr"), !! as.name("Position_physical"))
  }
  if (type == "map") {
    map_frame <- MF.obj$input[[map]]
    df <- inner_join(MF.obj$ref_map, map_frame, by = "Marker") %>%
      arrange(!! as.name("Chr"), !! as.name("Position_physical"))
  }
  MF.obj$config$removed_markers <- NULL
  MF.obj$pspline$recombination_rate <- NULL
  MF.obj$pspline$gam_models <- list()
  MF.obj$pspline$gam_fitted <- NULL
  MF.obj$pspline$gam_validate <- NULL

  for (i in seq_along(chromosomes)) {
    df_i <- filter(df, (!!quote(LG)) == chromosomes[i] & (!!quote(Chr)) == chromosomes[i])
    ## Through error if dataset is small
    if (nrow(df_i) <= 5) {
      stop(paste0("Not enough data points to fit Pspline for chromosome ", chromosomes[i]))
    }
    k <- (nrow(df_i) / 2) %>%
      trunc()
    if (k > 25) {
      k <- 25
    }
    set.seed(987428)
    MF.obj$pspline$gam_models[[i]] <- gam(
      Position ~ s(Position_physical, bs = "tp", fx = FALSE, k = k),
      data = df_i, method = "GCV.Cp"
    )
    df_i$Predictions <- predict(MF.obj$pspline$gam_models[[i]], df_i)
    markers_toremove <- which(residuals(MF.obj$pspline$gam_models[[i]], type = "scaled.pearson") %>%
      abs() > z)
    while (length(markers_toremove) != 0) {
      df_i <- df_i[-markers_toremove, ]
      MF.obj$config$removed_markers <- rbind(MF.obj$config$removed_markers, df_i[markers_toremove, ])
      MF.obj$pspline$gam_models[[i]] <- NULL
      MF.obj$pspline$gam_models[[i]] <- gam(
        Position ~ s(Position_physical, bs = "tp", fx = FALSE, k = k),
        data = df_i, method = "GCV.Cp"
      )
      df_i$Predictions <- predict(MF.obj$pspline$gam_models[[i]], df_i)
      markers_toremove <- NULL
    }
    MF.obj$pspline$gam_fitted <- rbind(MF.obj$pspline$gam_fitted, df_i)
    ## Fit in dense interval for later plotting
    df_pred <- data.frame(Marker = chromosomes[i],
                          LG = chromosomes[i],
      Position = seq(0,trunc(max(df_i$Position)), length.out = 250),
                          Position_physical = seq(0, trunc(max(df_i$Position_physical)),length.out = 250))
    df_pred$Predictions <- predict(MF.obj$pspline$gam_models[[i]], df_pred)
    MF.obj$pspline$gam_validate <- rbind(MF.obj$pspline$gam_validate, df_pred)
  }
  for (i in seq_along(chromosomes)) {
    rr_df <- filter(MF.obj$ref_map, (!!quote(Chr)) == chromosomes[i])
    r_rate <- data.frame(Position_physical = seq(min(rr_df$Position_physical), max(rr_df$Position_physical), by = 0.01))
    preds <- predict(MF.obj$pspline$gam_models[[i]], r_rate)
    names(preds) <- NULL
    r_rate$cM_per_Mbp <- c(diff(preds, differences = 1) * 10, NA)
    r_rate$Chr <- chromosomes[i]
    MF.obj$pspline$recombination_rate <- rbind(MF.obj$pspline$recombination_rate, r_rate)
  }
  MF.obj$pspline$recombination_rate$cM_per_Mbp[MF.obj$pspline$recombination_rate$cM_per_Mbp < 0] <- 0
  MF.obj$pspline$recombination_rate <- MF.obj$pspline$recombination_rate %>%
    na.omit()
  names(MF.obj$pspline$gam_models) <- chromosomes
  return(MF.obj)
}
