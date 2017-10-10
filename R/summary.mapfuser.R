# summary.mapfuser
#' internal position check
#'
#' @param object mapfuser object
#' @param ... ingored in function call
#' @return length and number of mapped marker for each input genetic map
#' @author Dennis van Muijen
#' @export

summary.mapfuser <- function(object, ...) {
  res <- lapply(object$raw_data, function(x) x %>%
    group_by("LG") %>%
    summarize(length = max(Position), n.markers = n()))
  print(res)
}
