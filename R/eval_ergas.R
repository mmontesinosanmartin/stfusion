#' Calculates the ERGAS index
#' 
#' Quantifies the error between the actual and predicted image
#' 
#' This is a frequently used evaluation metric for data fusion techniques.
#' ERGAS is the acronym for Erreur Relative Globale Adimensionnelle de Synthese.
#' A lower value means a lower error and therefore a better prediction. Both
#' \code{x} and \code{y}, should have the same extent and number of bands.
#' 
#' @param x the actual image as a \code{RasterStack}
#' @param y the predicted image as a \code{RasterStack}
#' 
#' @return the ERGAS value
#' 
#' @references Wald, L. (2002). Data fusion: definitions and architectures:
#' fusion of images of different spatial resolutions. Presses des MINES.
#' 
#' @example 
#' library(raster)
#' x1 <- x2 <- raster(nrow = 10, ncol = 10)
#' x1[] <- rnorm(100,100)
#' x2[] <- rnorm(100,100)
#' x <- stack(x1,x2)
#' y <- x
#' eval_ergas(x,y)
#'  
eval_ergas <- function(x,y){
  
  # info
  x.res <- res(x)
  y.res <- res(y)
  x.mat <- as.matrix(x)
  y.mat <- as.matrix(y)
  n.lyr <- ncol(x.mat)
  
  # computation
  res.r <- mean(x.res/y.res, na.rm = TRUE)
  avg.l <- apply(x.mat, 2, mean, na.rm = TRUE)
  rms.l <- sapply(1:n.lyr, function(i,x,y){
    .rmse(x[,i],y[,i])
  }, x = x.mat, y = y.mat)
  100 * res.r * sqrt(mean((rms.l/avg.l)^2, na.rm = TRUE))
  
}

.rmse <- function(x,y){
  sqrt(mean((x-y)^2, na.rm = TRUE))
}
