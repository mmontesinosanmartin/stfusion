#' Calculates the SSIM index
#' 
#' Quantifies the error between the actual and predicted image based on their
#' visual similarity
#' 
#' The SSIM acronym stands for Structural Similarity (Wang et al., 2004). The
#' index has three components, which are the structure, luminance and contrast.
#' The closer the SSIM is to 1, the more similar are the \code{x} and \code{y}.
#' Inner constants are kept as Wang et al., 2004 when \code{method = "default}.
#' If \code{method="uiqi"}, the function computes a special case of SSIM called
#' Universal Image Quality Index (UIQUI) (Wang and Bovik, 2002) with inner
#' constants equal to zero.
#' 
#' @param x the actual image as a \code{RasterStack}
#' @param y the predicted image as a \code{RasterStack}
#' @param scale the dynamic range of the image (default, \code{c(0,255)})
#' @param method 
#' 
#' @return a \code{list} with the SSIM value and the structural (\code{c}),
#' luminance (\code{l}), and contrast (\code{c}) evaluation components
#' 
#' @references
#' Wang, Z., & Bovik, A. C. (2002).
#' A universal image quality index.
#' IEEE signal processing letters, 9(3), 81-84.
#' 
#' Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004).
#' Image quality assessment: from error visibility to structural similarity.
#' IEEE transactions on image processing, 13(4), 600-612.
#' 
#' @example 
#' library(raster)
#' x1 <- x2 <- raster(nrow = 10, ncol = 10)
#' x1[] <- rnorm(100,100)
#' x2[] <- rnorm(100,100)
#' x <- stack(x1,x2)
#' y <- x
#' eval_ssim(x,y,scale=c(0,255))
#' 
eval_ssim <- function(x,y,scale=c(0,1), method = "default"){

  # stats
  x.mu <- cellStats(x, mean, na.rm = TRUE)
  y.mu <- cellStats(y, mean, na.rm = TRUE)
  x.sg <- cellStats(x, sd, na.rm = TRUE)
  y.sg <- cellStats(y, sd, na.rm = TRUE)
  xy.c <- diag(cov(as.matrix(x),as.matrix(y), use = "complete.obs" ))
  
  # constants
  if(method == "default"){
    c1 <- (0.01 * scale[2])^2
    c2 <- (0.03 * scale[2])^2
    c3 <- c2/2
  } else if(method == "uiqi") {
    c1 <- c2 <- c3 <- 0
  }

  # luminance
  lumi <- (2 * x.mu * y.mu + c1)/(x.mu^2 + y.mu^2 + c1)
    
  # contrast
  cont <- (2 * x.sg * y.sg + c2)/(x.sg^2 + y.sg^2 + c2)
  
  # structure
  stru <- (xy.c + c3)/((x.sg * y.sg) + c3)
  
  return(list(ssim=mean(stru * lumi * cont),
              s = mean(stru),
              l = mean(lumi),
              c = mean(cont)))
}
