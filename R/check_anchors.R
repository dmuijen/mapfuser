# check_anchors
#' Replace sub-linkage groups (e.g. 1.1, 1.2) with truncated linkage group number and split per chromosome for easy integration
#' @param MF.obj mapfuser object
#' @param anchors Number of uniquely positioned anchor markers between genetic maps
#' @return mapfuser object with QC$maps list item containing input maps splitted to unique linkage groups

check_anchors <- function(MF.obj = MF.obj, anchors = anchors) {
  anchor.markers <- list()
  if (anchors < 3) {
    stop("Minimum of three anchor markers required")
  }
  remove_list <- list()
  MF.obj$QC$networks <- NULL
  chr <- unique(bind_rows(MF.obj$QC$maps)$LG)
  for (i in seq_along(chr)) {
    idx <- which(unlist(lapply(MF.obj$QC$maps, function(x){
      unique(x[,"LG"] == chr[i])
    })))
    lg_maps <- MF.obj$QC$maps[idx] 
    lg_id <- lg_maps %>%
      names() %>%
      unique()
    anchor.markers_i <- matrix(
      vector(), length(lg_maps), length(lg_maps),
      dimnames = list(lg_id, lg_id)
    )
    for (l in seq_along(lg_maps)) {
      edge1 <- lg_maps[[l]]
      for (k in seq_along(lg_maps)) {
        anchor_df <- inner_join(lg_maps[[l]], lg_maps[[k]], by = "Marker", suffix = c("_map1", "_map2")) %>%
          group_by(!!as.name("Marker")) %>%
          slice(1)
        if (nrow(anchor_df) == 0) {
          anchor.markers_i[k, l] <- 0
        }
        if (nrow(anchor_df) != 0) {
          if (anchor_df$Position_map1 %>%
            unique() %>%
            length() < anchors | anchor_df$Position_map2 %>%
            unique() %>%
            length() < anchors) {
            anchor.markers_i[k, l] <- 0
          }
          if (nrow(anchor_df) >= 3) {
            condition1 <- any(sort(anchor_df$Position_map1 %>%
              diff() %>%
              abs(), decreasing = TRUE)[1:2] < 1)
            condition2 <- any(sort(anchor_df$Position_map2 %>%
              diff() %>%
              abs(), decreasing = TRUE)[1:2] < 1)
            if (any(condition1, condition2) == TRUE) {
              anchor.markers_i[k, l] <- 0
            } else {
              anchor.markers_i[k, l] <- nrow(anchor_df)
            }
          }
        }
      }
    }
    diag(anchor.markers_i) <- 0
    anchor.markers_i[anchor.markers_i < anchors] <- 0 
    MF.obj$QC$networks[[i]] <- graph.adjacency(anchor.markers_i %>%
      as.matrix(), mode = c("undirected"), weighted = TRUE)
    MF.obj$QC$networks[[i]] <- set_vertex_attr(MF.obj$QC$networks[[i]], name = "chromosome", value = chr[i])
    E(MF.obj$QC$networks[[i]])$weight <- 1 / E(MF.obj$QC$networks[[i]])$weight
    insuf_anchor <- NULL
    if (any(apply(anchor.markers_i, 1, max) < anchors)) {
      insuf_anchor <- which(apply(anchor.markers_i, 1, max) < anchors)
      if (length(insuf_anchor) > 0) {
        anchor.markers_i <- anchor.markers_i[-insuf_anchor, -insuf_anchor]
        for (z in seq_along(insuf_anchor)) {
          cat(paste0("Insufficient anchors..removed ", names(insuf_anchor)[z]), "\n")
        }
      }
    }
    remove_list[[i]] <- names(insuf_anchor)
    anchor.markers[[i]] <- anchor.markers_i
  }
  names(MF.obj$QC$networks) <- chr
  ### Remove maps with insufficient linkage to at least one other map
  remove_list <- unlist(remove_list)
  if (!is.null(remove_list)) {
    MF.obj$config$removed_LGIDs <- c(MF.obj$config$removed_LGIDs, remove_list)
    MF.obj$QC$maps <- MF.obj$QC$maps[!names(MF.obj$QC$maps) %in% remove_list]
  }

  ### Return QC corrected networks for MST
  MF.obj$QC$qcnetworks <- MF.obj$QC$networks
  if (!is.null(remove_list)) {
    for (i in seq_along(MF.obj$QC$qcnetworks)) {
      rm_id <- remove_list[remove_list %in% V(MF.obj$QC$qcnetworks[[i]])$name]
      MF.obj$QC$qcnetworks[[i]] <- delete_vertices(MF.obj$QC$qcnetworks[[i]], rm_id)
    }
  }
  return(MF.obj)
}
