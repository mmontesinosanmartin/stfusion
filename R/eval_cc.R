#' Calculates the Pearson's Correlation Coefficient (CC)
#' 
#' Quantifies the average pair-wise correlation between the bands of two images
#' 
#' \code{x} and \code{y} are assumed to have the same resolution and number of
#' bands.
#' 
#' @param x the high-resolution actual image as a \code{RasterStack}
#' @param y the high-resolution predicted image as a \code{RasterStack}
#' @param by.band whether to compute the statistic by layer
#' @param n.smpls number of samples to compute the statistic (default, \code{NULL})
#' 
#' @example 
#' x <- raster(nrow=10,ncol=10)
#' x[] <- rnorm(100)
#' y <- raster(nrow=10,ncol=10)
#' y[] <- rnorm(100)
#' eval_cc(x,y)
#' 
eval_cc <- function(x, y, by.band = FALSE, n.smpl = NULL){
  
  # check projection, extent and resolution
  x.siz <- prod(res(x))
  y.siz <- prod(res(projectRaster(raster(y), crs = projection(x))))
  mxr <- x.siz < y.siz
  ref <- x; if(mxr) ref <- y
  slv <- y; if(mxr) slv <- x
  slv <- as(st_warp(st_as_stars(slv),
                 st_as_stars(ref),
                 method = "average",
                 use_gdal = TRUE), "Raster")
  npx <- prod(dim(ref)[1:2])
  
  # sampling?
  if(is.null(n.smpl)){
    smpl <- 1:npx
  } else {
    smpl <- sample(1:npx, n.smpl)
  }
  
  # Compute correlation
  cors <- diag(cor(as.matrix(ref[])[smpl,],
                   as.matrix(slv[])[smpl,],
                   use = "complete.obs"))
  
  # format of the output
  if(by.band){
    out <- cors
  } else {
    out <- mean(cors, na.rm = TRUE)
  }
  return(out)
}
