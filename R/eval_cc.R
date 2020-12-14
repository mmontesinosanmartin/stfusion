#' Calculates the Pearson's Correlation Coefficient (CC)
#' 
#' Quantifies the average pair-wise correlation between the bands of two images
#' 
#' \code{x} and \code{y} are assumed to have the same resolution and number of
#' bands.
#' 
#' @param x the high-resolution actual image as a \code{RasterStack}
#' @param y the high-resolution predicted image as a \code{RasterStack}
#' 
#' @example 
#' x <- raster(nrow=10,ncol=10)
#' x[] <- rnorm(100)
#' y <- raster(nrow=10,ncol=10)
#' y[] <- rnorm(100)
#' eval_cc(x,y)
#' 
eval_cc <- function(x, y, by.band = FALSE){
  cors <- diag(cor(x[], y[],use = "complete.obs"))
  if(by.band){
    out <- cors
  } else {
    out <- mean(cors, na.rm = TRUE)
  }
  return(out)
}
