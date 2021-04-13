#' Computes the Root Mean Square Error (RMSE)
#' 
#' Quantifies the RMSE between the bands of two images
#' 
#' \code{x} and \code{y} are assumed to have the same projection, spatial
#' resolution, and number of bands.
#' 
#' @param x the high-resolution actual image as a \code{RasterStack}
#' @param y the high-resolution predicted image as a \code{RasterStack}
#' @param by.band if \code{FALSE}(default), gives the average RMSE of all bands
#' 
#' @example 
#' x <- raster(nrow=10,ncol=10)
#' x[] <- rnorm(100)
#' y <- x
#' eval_rmse(x,y)
#' 
eval_rmse <- function(x, y, by.band = FALSE, rm.na = TRUE){
  x.mat <- as.matrix(x[])
  y.mat <- as.matrix(y[])
  compl <- 1:nrow(y.mat)
  if(rm.na) compl <- complete.cases(cbind(x.mat, y.mat))
  out <- cpp_rmse(x.mat[compl,], y.mat[compl,], byband = by.band)
  return(out)
  # rmse.bnd <- sqrt(apply(as.matrix(((x-y)^2)[]),2, mean, na.rm = TRUE))
  # if(by.band){
  #   rmse.bnd
  # } else {
  #   mean(rmse.bnd, na.rm = TRUE)
  # }
}
