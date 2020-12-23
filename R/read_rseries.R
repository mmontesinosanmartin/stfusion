#' @title Reads several multi-spectral images
#' 
#' @description reads several multi-spectral images from their file paths
#' and returns a list of \code{RasterStack}a. A wrapper function of 
#' \code{stack()}.
#' 
#' @details This is a function to handle series of multi-spectral
#' images using \code{raster}.
#' 
#' @return a list of \code{RasterStacks} named by dates
#' 
#' @param filenames a vector of file paths
#' 
read_rseries <- function(filenames){
  
  dates <- get_dates(basename(filenames), "%Y%j")
  imgs <- lapply(filenames, raster::stack)
  names(imgs) <- dates
  
  return(imgs)
}
