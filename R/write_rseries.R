#' Writes a series of images
#' 
#' @param x a \code{List} with raster time series
#' @param filenames a \code{character} vector of file names
#' @param overwrite whether to overwrite existing files. Default \code{TRUE} 
#' 
#' @return nothing is returned

write_rseries <- function(x, filenames, overwrite = TRUE){
  
  nimgs <- length(x)
  
  for(i in 1:nimgs){
    
    writeRaster(x[[i]], filename = filenames[i], overwrite = overwrite)
    
  }
}
