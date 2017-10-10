# LPmerge_par
#' Wrapper for multicore and multichromosome merging of maps using the LPmerge
#'
#' The LPmerge algorithm is called for all linkage groups with an option for multicore processing. The mapfuser package further adds automatic selection of max.interval giving the lowest RMSE
#' @param MF.obj A mapfuser object genetics maps loaded and optionally a reference map
#' @param max.interval A whole number specifying the maximum interval size between bins to include in the objective function. An array of numbers can be passed to test different values (one consensus map is produced for each value in the array).
#' @param weights Optional vector of length T containing the weights for each map in the objective function (see details). If not passed, the maps are given equal weight.
#' @param max.int_sel Either automatically select the max.interval with the lowest RMSE for each linkage group or specify manually a vector of values of max.interval to select for the output
#' @param n.cores number of cores
#' @return Three items added to the mapfuser object under the results slot. 1) The consensus map at the selected max.interval, either manually specified or automatically select
#' 2) A list of calculated Root Mean Square Error (RMSE) per chromosome with RMSE for each linkage group ID compared to the consensus map. THe mean and standard deviation over the linkage group IDs is also provided.
#' 3) A list with length equal to the length of the max.interval parameter. Each entry in the list is a data frame containing the consensus map and the component linkage maps.
#' @import LPmerge doParallel foreach parallel
#' @references Endelman, JB, and C Plomion. 2014. LPmerge: An R package for merging genetic maps by linear programming. Bioinformatics 30:1623-1624.
#' @author Dennis van Muijen
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "-1", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- map_qc(MF.obj)
#' MF.obj$config$chr <- "1"
#' MF.obj <- LPmerge_par(MF.obj = MF.obj, n.cores = 1,
#' max.interval = 1, max.int_sel = "auto", weights = NULL)
#' }
#' \dontrun{
#' MF.obj <- LPmerge_par(MF.obj = MF.obj, n.cores = 2,
#' max.interval = 1, max.int_sel = "auto", weights = NULL)
#' # Plot result
#' plot(MF.obj, which = "single_map", maps = "consensus")
#' # Access RMSE table
#' MF.obj$result$RMSE
#' }
#' @export

LPmerge_par <- function(MF.obj, n.cores = 2, max.interval = 1:3, max.int_sel = "auto", weights = NULL) {
  if (class(MF.obj) != "mapfuser") {
    stop("Object should be of class mapfuser")
  }
  if (is.null(MF.obj$QC)) {
    stop("Run Quality Control First")
  }
  if (MF.obj$config$LPmerge_ready == FALSE) {
    stop("Run map_qc succesfully first")
  }
  i <- NULL
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
  chr <- MF.obj$config$chr
  lp_res <- foreach(i = seq_along(chr), .packages = c("dplyr", "Matrix", "Rglpk", "slam", "LPmerge"), .inorder = TRUE) %dopar% {
    lg <- lapply(MF.obj$QC$maps, function(x) filter(x, (!!quote(LG)) == chr[i]))
    lg <- lg[sapply(lg, function(x) dim(x)[1]) > 1] 
    lg_map <- lapply(lg, "[", TRUE, c(2, 4)) 
    weights <- lapply(lg, function(x) x$Weight %>%
      unique()) %>%
      unlist()
    if (length(lg_map) > 1) {
      consensusmap <- LPmerge(Maps = lg_map, max.interval = max.interval, weights = weights)
    } else {
      consensusmap <- "Nothing to merge.."
    }
  }
  mapnames <- NULL
  for (k in seq_along(chr)) {
    if (!is.null(lp_res[[k]])) {
      df <- lapply(MF.obj$QC$maps, function(x) filter(x, (!!quote(LG)) == chr[k]))
      df <- df[sapply(df, function(x) dim(x)[1]) > 1]
      mapnames <- lapply(df, function(x) {
        x$LGID %>%
          unique()
      })
      colnames(lp_res[[k]][[1]]) <- c("marker", "consensus", do.call(rbind, mapnames)[, 1]) %>%
        as.character()
    }
  }
  stopCluster(cl)
  map_RMSE <- calc_RMSE(lp_res, max.interval = max.interval)
  consensus_map <- select_result(lp_res = lp_res, map_RMSE = map_RMSE, max.int_sel = max.int_sel, chr = chr)
  result_list <- vector("list", length = 3)
  if (max.int_sel[1] == "auto") {
    max.int.index <- lapply(map_RMSE, function(x) {
      which.min(tail(x, 2)[1, ][, -1]) %>%
        names()
    }) %>%
      unlist()
    MF.obj$config$max.interval <- max.int.index
  } else {
    MF.obj$config$max.interval <- max.int_sel
  }
  MF.obj$result$consensus_map <- consensus_map
  names(MF.obj$result$consensus_map) <- c("Marker", "LG", "Position")
  MF.obj$result$consensus_map$LG <- MF.obj$result$consensus_map$LG %>%
    as.factor()
  MF.obj$result$RMSE <- map_RMSE
  MF.obj$result$LPmerge <- lp_res
  MF.obj$config$selected_max.int <- attributes(MF.obj$result$consensus_map)$max.interval
  return(MF.obj)
}

# calc_RMSE
#' Internal function to calculate the Root Mean Square Error to select the maximum interval size.
#'
#' The LPmerge algorithm is called for all linkage groups with an option for multicore processing. mapfuser further adds automatic selection of the lowest max
#' @param lp_res Result of the mapfuser foreach call or alternatively the result of a call to LPmerge
#' @param max.interval A whole number specifying the maximum interval size between bins to included in the objective function for LPmerge
#' @return A list of calculated Root Mean Square Error (RMSE) per chromosome with RMSE for each linkage group ID compared to the consensus map. The mean and standard deviation over the linkage group IDs is also provided.
#' @author Dennis van Muijen

calc_RMSE <- function(lp_res, max.interval) {
  lapply(lp_res, function(x) {
    RMSE <- lapply(x, function(y) {
      apply(y[, -1], 2, function(z) {
        sqrt(mean((y$consensus - z) ^ 2, na.rm = TRUE))
      })
    })
    df <- data.frame(map = names(RMSE[[1]]))
    row.names(df) <- NULL
    for (z in seq_along(RMSE)) {
      df <- cbind(df, max.int = RMSE[[z]])
    }
    colnames(df) <- c("map", paste0(rep("max.interval", length(RMSE)), max.interval))
    df <- df[-1, ]
    df$map <- df$map %>%
      as.character()
    if (ncol(df) > 2) {
      df <- rbind(df, c("mean", apply(df[, -1], 2, mean)))
      df <- rbind(df, c("sd", apply(df[, -1], 2, sd))) %>%
        as.data.frame()
      row.names(df) <- NULL
      df[, -1] <- apply(df[-1], 2, as.numeric)
    } else {
      df <- rbind(df, c("mean", mean(df[, -1])))
      df <- rbind(df, c("sd", sd(df[, -1]))) %>%
        as.data.frame()
      row.names(df) <- NULL
      df[, -1] <- as.numeric(df[, -1])
    }
    return(df)
  })
}

# select_result
#' Select result from mapfuser
#'
#' Internal function for mapfuser -  Flatten, simplify, and combine to one dataframe the consensus map of all linkage groups for easy inspection and exporting
#' @param lp_res Result of the mapfuser foreach call or alternatively the result of a call to LPmerge
#' @param max.int_sel Either automatic selection of the maximum interval size that minimized the mean Root Square Mean Error between the consensus map and invidual linkage maps
#' @param chr chromosomes to perform analysis for
#' @param map_RMSE Output of a call to calc_RMSE
#' @return The consensus map at the selected maximum interval size K in a convenient format
#' @author Dennis van Muijen

select_result <- function(lp_res, max.int_sel = "auto", map_RMSE, chr = chr) {
  if (max.int_sel[1] == "auto") {
    cat("Selected max.interval with lowest RMSE..", fill = TRUE)
    max.int.index <- lapply(map_RMSE, function(x) {
      which.min(tail(x, 2)[1, ][, -1])
    })
    selected_maps <- list()
    for (u in seq_along(lp_res)) {
      to_select <- max.int.index[[u]]
      selected_maps[[u]] <- lp_res[[u]][[to_select]][,1:2]
      selected_maps[[u]] <- data.frame(chromosome = chr[u], selected_maps[[u]], stringsAsFactors = F)
    }
    consensus_map <- bind_rows(selected_maps)[, c(2, 1, 3)]
    attributes(consensus_map)$max.interval <- max.int.index %>%
      unlist() %>%
      as.vector()
    return(consensus_map)
  }
  if (max.int_sel[1] != "auto") {
    if (length(max.int_sel) != length(lp_res)) {
      stop("wrong max.interval values specified")
    }
    selected_maps <- list()
    for (v in seq_along(lp_res)) {
      to_select <- max.int_sel[v]
      selected_maps[[]] <- lp_res[[v]][[to_select]]
      selected_maps[[v]] <- cbind(chromosome = chr[v], selected_maps[[v]])
    }
    consensus_map <- bind_rows(selected_maps)[, c(2, 1, 3)]
    attributes(consensus_map)$max.interval <- max.int_sel
    return(consensus_map)
  }
}