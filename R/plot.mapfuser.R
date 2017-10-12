# plot.mapfuser
#' Visualise mapfuser object data
#'
#' @param x A mapfuser object with a reference map loaded and either a consensus map created with mapfuser or a genetic loaded with read_maps
#' @param which Dataset to visualize:
#' 1) which = "mapnetwork", "mst". Plots the connection in terms of overlapping markers between maps as mapnetwork and the corresponding minimum spanning tree using the inverse of the number of anchors markers as edge weight
#' 2) which = "compare_maps". Plots a simple scatterplot between two genetic maps,
#' 3) which = "single_map". Plot a single map per chromosome in mapchart style
#' 4) which = "marey_map". Plot the relationship between centi Morgan positions and physical genome position along with the fitted Penalized spline for each chromosomes
#' 5) which = "recombination_rate". Plot the recombination rate in centi Morgan per mega base pair along the physical genome.
#' @param chr The chromosomes to plot, mapnetwork and mst is single chromosome only
#' @param maps Name or names of the maps to plot in case which = "mareymap", "single_map", or "compare_maps"
#' @param ... ingored in function call
#' @return A plot.igraph object in the case of "mst" or "mapnetwork" or an interactive ggplotly object otherwise
#' @author Dennis van Muijen
#' @examples
#' \dontshow{
#' fpath <- system.file("extdata", package="mapfuser")
#' maps <- list.files(fpath, pattern = "-1", full.names = TRUE)
#' MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
#' MF.obj <- map_qc(MF.obj)
#' plot(x = MF.obj, which = "mapnetwork", chr = 1)
#' }
#' \dontrun{
#' plot(x = MF.obj, which = "mapnetwork", chr = 1)
#' plot(x = MF.obj, which = "mst", chr = 1)
#' plot(x = MF.obj, which = "single_map", maps = "consensus")
#' plot(x = MF.obj, which = "compare_maps", maps = c("Col-0_Cvi-0.csv","Col-0_Sha.csv"), chr = 1:3)
#' plot(x = MF.obj, which = "mareymap", maps = "Col-0_Bur-0.csv", chr = 1:5)
#' plot(x = MF.obj, which = "mareymap", maps = "consensus", chr = 1:5)
#' plot(x = MF.obj, which = "recombination_rate", chr = 1:5)
#' }
#' @export

plot.mapfuser <- function(x, which = c("mapnetwork", "mst", "compare_maps", "single_map", "mareymap", "recombination_rate"),
                          chr = NULL, maps = NULL, ...) {
  if (class(x) != "mapfuser") {
    stop("Object should be of class mapfuser")
  }
  if (which == "mapnetwork") {
    chrom <- x$config$chr[x$config$chr == chr] %>%
      as.character()
    G <- x$QC$networks[[chrom]]
    p <- visIgraph(G)
    return(p)
  }
  if (which == "mst") {
    chrom <- x$config$chr[x$config$chr == chr] %>%
      as.character()
    G <- x$QC$mst[[chrom]]
    p <- visIgraph(G)
    return(p)
  }
  if (which == "compare_maps") {
    if (is.null(maps) | length(maps) != 2) {
      stop("Specify chr and/or two maps to compare")
    }
    df <- inner_join(x$input[[maps[1]]], x$input[[maps[2]]], by = "Marker")
    names(df) <- gsub("[.]", "_", names(df))
    if (!is.null(chr)) {
      df <- filter(df, (!!quote(LG_x)) %in% chr)
    }
    p <- ggplot(aes_string(x = "Position_x", y = "Position_y"), data = df)
    p <- p + geom_point() + facet_wrap("LG_x", scales = "free") + theme_bw() +
      xlab(paste0(maps[1], " cM")) + ylab(paste0(maps[2], " cM"))
    return(ggplotly(p))
  }
  if (which == "mareymap") {
    if (maps != "consensus") {
      if (is.null(maps) | length(maps) != 1) {
        stop("Specify maps argument")
      }
      if (is.null(x$ref_map)) {
        stop("No reference map loaded")
      }
      df <- inner_join(x$ref_map, x$input[[maps[1]]], by = "Marker")
      names(df) <- gsub("[.]", "_", names(df))
      if (is.null(chr)) {
        chromosomes <- x$config$chromosomes
      } else {
        chromosomes <- chr
      }
      p <- ggplot(aes_string(x = "Position_physical", y = "Position", data_id = "Marker"), data = filter(x$pspline$gam_fitted, (!!quote(LG)) %in% chromosomes))
      p <- p + geom_point() + facet_wrap("LG", scales = "free") + theme_bw() +
        xlab(x$config$ref_name) + ylab("Consensus map (cM)") +
        geom_line(aes_string(x = "Position_physical", y = "Predictions", group = 1), colour = "red", lwd = 0.4, data = filter(x$pspline$gam_validate, (!!quote(LG)) %in% chromosomes)) +
        guides(colour = FALSE)
      return(ggplotly(p))
    }
    if (maps == "consensus") {
      if (is.null(maps) | length(maps) != 1) {
        stop("Specify maps argument")
      }
      if (is.null(x$pspline)) {
        stop("Run genphys_fit first")
      }
      if (is.null(chr)) {
        chromosomes <- x$config$chromosomes
      } else {
        chromosomes <- chr
      }
      p <- ggplot(aes_string(x = "Position_physical", y = "Position", data_id = "Marker"), data = filter(x$pspline$gam_fitted, (!!quote(LG)) %in% chromosomes))
      p <- p + geom_point() + facet_wrap("LG", scales = "free") + theme_bw() +
        xlab(x$config$ref_name) + ylab("Consensus map (cM)") +
        geom_line(aes_string(x = "Position_physical", y = "Predictions", group = 1), colour = "red", lwd = 0.4) +
        guides(colour = FALSE)
      return(ggplotly(p))
    }
  }
  if (which == "recombination_rate") {
    if (is.null(x$pspline$recombination_rate)) {
      stop("Run genphys_fit first")
    }
    if (is.null(chr)) {
      chromosomes <- x$config$chromosomes
    } else {
      chromosomes <- chr
    }
    p <- ggplot(aes_string(x = "Position_physical", y = "cM_per_Mbp"), data = filter(x$pspline$recombination_rate, (!!quote(Chr)) %in% chromosomes))
    p <- p + geom_line(group = 1, lwd = 1.4) + facet_wrap("Chr", scales = "free_x") + theme_bw()
    return(ggplotly(p))
  }
  if (which == "single_map") {
    if (is.null(maps)) {
      stop("Specify maps argument")
    }
    if (maps == "consensus") {
      ggdata <- x$result$consensus_map
    } else {
      ggdata <- x$input[[maps[1]]]
    }
    p <- ggplot(aes_string(x = "LG", y = "Position", data_id = "Marker"), data = ggdata) +
      geom_path(data = ggdata, aes_string(x = "LG", y = "Position", group = "LG"), size = 4.5, lineend = "round", colour = "grey") +
      theme_bw() + ylab("cM") +
      xlab("Linkage Group") + scale_y_reverse() +
      geom_point(size = 3.5, col = "black", shape = 95)
    return(ggplotly(p = p, tooltip = c("x", "data_id")))
  }
}
