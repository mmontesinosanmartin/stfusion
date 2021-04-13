#' Calculates the QNR index
#' 
#' Quantifies the error between the actual and predicted image in terms
#' of spectral and spatial quality measurements
#' 
#' QNR stands for Quality with No Reference. Its current implementation is an 
#' adaptation from Alparone et al., (2008). The QNR is a measurement of
#' distortion where values closer to 1 indicate better predictions. The
#' function assumes that \code{x.h}, \code{x.l}, and \code{y} have the
#' same number of pixels and projection, i.e. they low resolution image
#' has been re-projected and re-sampled.
#' 
#' @param x.h the high-resolution actual image as a \code{RasterStack}
#' @param x.l the low-resolution actual image as a \code{RasterStack}
#' @param y the predicted image as a \code{RasterStack}
#' @param scale the dynamic range of the image (default, \code{c(0,255)})
#' 
#' @return the QNR value
#' 
#' @references Alparone, L., Aiazzi, B., Baronti, S., Garzelli, A., Nencini, F., & Selva, M. (2008).
#' Multispectral and panchromatic data fusion assessment without reference.
#' Photogrammetric Engineering & Remote Sensing, 74(2), 193-200.
#' 
#' @example 
#' library(raster)
#' x1 <- x2 <- raster(nrow = 9, ncol = 9)
#' x1[] <- rnorm(81,100)
#' x2[] <- rnorm(81,100)
#' x.h <- stack(x1,x2)
#' x.l <- aggregate(x.h, fact=3)
#' y <- x.h
#' eval_aqnr(x.h,x.l,y)
#' 
eval_qnr <- function(x.h,x.l,y,scale=c(0,255)){
  
  # low-res versions
  n.lyrs <- nlayers(x.h)
  
  # standardization
  x.l <- projectRaster(x.l, y, method = "ngb")
  # x.h <- projectRaster(x.h, y, method = "ngb")
  
  # spectral
  spec.diff <- rep(NA, n.lyrs * (n.lyrs - 1))
  counter <- 0
  for(i in 1:n.lyrs){
    for(j in (1:n.lyrs)[-i]){
      counter <- counter + 1
      spec.diff[counter] <- abs(eval_ssim(x.l[[i]], x.l[[j]], scale = scale)$ssim -
                                eval_ssim(y[[i]], y[[j]], scale = scale)$ssim)
    }
  }
  d.lambda <- sqrt(mean(spec.diff))
  
  # spatial
  spat.diff <- rep(NA, n.lyrs)
  for(i in 1:n.lyrs){
    spat.diff[i] <- abs((eval_ssim(x.l[[i]], x.h[[i]], scale=scale)$ssim -
                         eval_ssim(x.h[[i]], y[[i]], scale=scale)$ssim))
  }
  d.spatil <- sqrt(mean(spat.diff))
  
  return((1 - d.spatil) * (1 - d.lambda))
  
}
