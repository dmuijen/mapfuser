# map_orient
#' Corrects orientation of input genetic maps
#' @param MF.obj A mapfuser object
#' @return A mapfuser object with corrected map orientation
#' @author Dennis van Muijen
#' @import mgcv igraph
#' @importFrom dplyr select filter group_by distinct arrange mutate summarize n inner_join bind_rows slice

map_orient <- function(MF.obj) {
  LG <- NULL 
  LGID <- NULL
  Chr <- NULL
  chr <- unique(bind_rows(MF.obj$QC$maps)$LG)
  edge.list <- list()
  for (i in seq_along(chr)) {
    edge.list[[i]] <-
      ends(MF.obj$QC$mst[[i]], es = E(MF.obj$QC$mst[[i]])) %>%
      as.data.frame(stringsAsFactors = FALSE)
  }
  ### Add map inversion if stepwise correlation is negative
  for (i in seq_along(chr)) {
    lg_maps <- lapply(MF.obj$QC$maps, filter, (!!quote(LG)) == chr[i])
    lg_maps <- lg_maps[lapply(lg_maps,nrow)>0]
    for (j in 1:nrow(edge.list[[i]])) {
      cor_index <- which(lapply(lg_maps, function(x)
        x$LGID %>%
          unique()) %>%
          unlist() %in% edge.list[[i]][j,])
      LGIDs <- (lapply(lg_maps, function(x)
        x$LGID %>%
          unique()) %>%
          unlist())[cor_index]
      maps_tojoin <-
        lg_maps[lapply(lg_maps, function(x)
          x$LGID %in% LGIDs %>%
            unique()) %>%
            unlist()]
      df <-
        inner_join(maps_tojoin[[1]], maps_tojoin[[2]], by = "Marker")
      edge.list[[i]][j, 3] <-
        cor(df$Position.x, df$Position.y, use = "pairwise.complete.obs")
    }
  }
  ### Extract path covering the mst
  connections <- vector("list", length(MF.obj$QC$mst))
  for (i in seq_along(MF.obj$QC$mst)) {
    full_path <- dfs(MF.obj$QC$mst[[i]], 1, unreachable = T)
    for (j in 1:(length(full_path$order) - 1)) {
      connections[[i]][[j]] <-
        all_simple_paths(
          MF.obj$QC$mst[[i]],
          from = full_path$order[j] %>%
            as.numeric(),
          to = full_path$order[j + 1] %>%
            as.numeric()
        )[[1]]
    }
    z_index <- NULL
    z <- connections[[i]] %>%
      unlist()
    for (k in 1:(length(z) - 1)) {
      if (z[k] == z[k + 1]) {
        z_index[k] <- k
      }
    }
    if (!is.null(z_index)) {
      connections[[i]] <- z[-z_index %>%
                              na.omit()]
    } else {
      connections[[i]] <- z
    }
  }
  full_map <- bind_rows(MF.obj$QC$maps)
  inverted_maps <- NULL
  for (i in seq_along(chr)) {
    for (j in 1:(length(connections[[i]]) - 1)) {
      mst_stepwise <- names(connections[[i]])[j:(j + 1)]
      lg_map <-
        filter(full_map, (!!quote(LGID)) %in% mst_stepwise)
      lg_map <- split(lg_map, lg_map$id)
      df <- inner_join(lg_map[[1]], lg_map[[2]], by = "Marker")
      cor_map <-
        cor(df$Position.x, df$Position.y, use = "pairwise.complete.obs")
      if (cor_map < 0) {
        to_invert <- names(connections[[i]][j + 1])
        full_map[which(full_map$LGID == to_invert),] <-
          map_flip(full_map[which(full_map$LGID == to_invert),])
        print(paste0("Inverted linkage group in map.. ", to_invert))
        inverted_maps <- c(inverted_maps, to_invert)
      }
    }
  }
  MF.obj$config$inverted_LGids <- inverted_maps
  if (is.null(MF.obj$ref_map)) {
    MF.obj$QC$maps <- split(full_map, full_map$LGID)
    return(MF.obj)
  }
  list_qcmaps <- split(full_map, full_map$LGID)
  ref_overlap <- vector("list", length(chr))
  for (i in seq_along(chr)) {
    cleaned_maps_lg <-
    lapply(list_qcmaps, filter, (!!quote(LG)) == chr[i])
    ref_overlap_i <- lapply(cleaned_maps_lg, function(x) {
      sum(x$Marker %in% MF.obj$ref_map$Marker)
    })
    ref_overlap[[i]] <- names(which.max(unlist(ref_overlap_i)))
  }
  # Check correlations
  sign <- NULL
  ## Take map with most anchors to reference map per chromosome
  ## Check if sufficient data
  for (i in seq_along(chr)) {
    df <-
      inner_join(
        MF.obj$ref_map,
        filter(
          full_map,
          (!!quote(LG)) == chr[i] &
            (!!quote(LGID)) == ref_overlap[[i]]
        ),
        by = "Marker"
      )
    if (nrow(df) == 0) {
      stop(paste0(
        "No overlap between reference map and linkage map for chromosome ",
        chr[i]
      ))
    }
    df <-
      arrange(filter(df, (!!quote(Chr)) == chr[i]),
              !!as.name("Position_physical"))
    if (nrow(distinct(df)) < 5) {
      cat(
        paste0(
          "Insufficient data to fit Pspline.. switching to correlation for chromosome "
        ),
        chr[i],
        "\n"
      )
      sign[i] <-
        cor(df$Position_physical, df$Position, use = "pairwise.complete.obs")
    }
    ### Overrule automatic model selection for small datasets..
    if (nrow(df) < 20) {
      k <- nrow(df) - 3
      k <- round(k, digits = 0)
    }
    if (nrow(df) >= 20) {
      k <- nrow(df) - 10
    }
    gam_fit <- gam(
      Position ~ s(
        Position_physical,
        bs = "tp",
        fx = FALSE,
        k = k
      ),
      data = df,
      method = "REML"
    )
    gam_pred <- predict(gam_fit, df)
    q <- diff(gam_pred) %>%
      mean()
    if (q > 0) {
      sign[i] <- 1
    }
    if (q <= 0) {
      sign[i] <- -1
    }
  }
  for (i in seq_along(chr)) {
    list_qcmaps <- lapply(list_qcmaps, function(x) {
      if (x$LG %>%
          unique() == chr[i] & sign[i] < 0) {
        x <- x %>%
          map_flip()
        return(x)
      } else {
        return(x)
      }
    })
  }
  inverted <- which(sign < 0)
  if (length(inverted) > 0) {
    for (i in seq_along(inverted)) {
      cat(paste0("Alligned linkage group ", chr[chr %in% inverted[i]], " to reference map"),
          sep = "\n")
    }
  }
  MF.obj$config$inverted_to_reference <- chr[inverted]
  MF.obj$QC$maps <- list_qcmaps
  return(MF.obj)
}