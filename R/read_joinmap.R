# read_joinmap
#' Internal function, read genetic maps in JoinMap format
#'
#' @param mapfile JoinMap formatted file to read
#' @return A dataframe with markers, linkage groups identifiers and genetic position on separate columns
#' @author Dennis van Muijen
#' @importFrom stringi stri_split_fixed
#' @importFrom tidyr separate
#' @importFrom dplyr bind_rows

read_joinmap <- function(mapfile){
  inFile <- read.delim(mapfile, comment.char = ';', header = FALSE, sep = "\t",
                       strip.white = TRUE, stringsAsFactors = FALSE)
  fix_spaces <- data.frame(mapvec = gsub("\\s+", " ",inFile[,1]))
  map_df <- separate(fix_spaces, !!quote(mapvec), sep = " ", into = c("Marker","Position"))
  idx <- grep("group",map_df$Marker)
  grp_idx <- stri_split_fixed(map_df$Position[idx]," ") %>% 
    unlist
  LGvec <- c()
  map_list <- vector("list", length(grp_idx))
  for(i in seq_along(grp_idx)){
    if(i < length(grp_idx)){
      map_list[[i]] <- map_df[idx[i]:(idx[(i+1)]-1),]
    }
    if(i == length(grp_idx)){
      map_list[[i]] <- map_df[idx[i]:nrow(map_df),]
    }
  }
  map_list_final <- lapply(map_list, function(x){
    x$LG <- x$Position[1]
    x <- x[-1,]
    return(x)
  })
  map_frame <- bind_rows(map_list_final)
  map_frame <- map_frame[,c(1,3,2)]
  return(map_frame)
}
