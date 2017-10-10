# map_split
#' Replace sub-linkage groups (e.g. 1.1, 1.2) with truncated linkage group number and split per chromosome for easy integration
#' @param MF.obj A mapfuser object
#' @return mapfuser object with QC$maps list item containing input maps splitted to unique linkage groups

map_split <- function(MF.obj) {
  '.' <- NULL
  MF.obj$QC$maps <- MF.obj$input %>%
    bind_rows(.id = "id") %>%
    mutate(LGID = paste0(.$id, "_", .$LG)) %>%
    arrange(!! as.name("id"),!! as.name("LG"), !! as.name("Position")) %>%
    mutate(LG = gsub("[.]1|[.]2|[.]3|[.]4|[.]5|[.]6|[.]7|[.]8|[.]9", "",!! as.name("LG"))) %>%
    split(.$LGID)
  return(MF.obj)
}