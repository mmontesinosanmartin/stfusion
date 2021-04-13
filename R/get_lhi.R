#' Calculates the Landscape Heterogeneity Index (LHI)
#' 
#' Quantifies the spatial hetergeneity of an image
#' 
#' The ``ddt'' stands for Direct Difference Threshold. This method computes
#' the difference between neighbouring pixels. A threshold determines whether
#' they are diffrent or not. According to Chen and Xu (2014), the ``ddt'' may
#' lead to biased estimates of spatial heterogeneity when pixel values are low.
#' 
#' As a result, Chen and Xu (2014) also propose the ``sp'' method, which stands
#' for Slope Projection. Instead of the differences, this method calculates the
#' ratio between neighboring pixels.
#' 
#' The LHI index ranges between 0 and 1, where 0 means homogeneus and 1
#' heterogeneous.
#' 
#' @param r a \class{RasterLayer}
#' @param method either "ddt" or "sp". See details
#' 
#' @return a number
#' 
#' @references Chen, B., & Xu, B. (2014). A novel method for measuring 
#' landscape heterogeneity changes. IEEE Geoscience and Remote Sensing
#' Letters, 12(3), 567-571.
#' 
#' @example 
#' library(raster)
#' r <- raster(nrow = 10, ncol = 10)
#' r[] <- rnorm(10,100)
#' get_lhi(r)
#' 
get_lhi <- function(r, method = "sp"){
  
  # to matrix
  r.mat <- as.matrix(r)
  n.col <- ncol(r.mat)
  n.row <- nrow(r.mat)
  
  # stats
  sd.img <- sd(r.mat, na.rm = T)
  mu.img <- mean(r.mat, na.rm = T)
  
  # ddt method
  if(method == "ddt"){

    # threshold
    nu.img <- 5
    th.img <- sd.img/nu.img
    
    # differences
    rv.img <- abs(apply(r.mat, 2, diff))
    rh.img <- abs(apply(r.mat, 1, diff))

  # sp method
  } else if(method == "sp"){
    
    # threshold
    th.img <- abs(sd.img/mu.img)
    
    # ratio
    rv.img <- abs(apply(r.mat, 2, function(x, n.row){x[2:n.row]/x[1:(n.row - 1)]},n.row))
    rh.img <- abs(apply(r.mat, 1, function(x, n.col){x[2:n.col]/x[1:(n.col - 1)]},n.col))

  }
  
  # apply threshold
  rv.img <- 1 * (rv.img > th.img)
  rh.img <- 1 * (rh.img > th.img)
  
  # final score
  (sum(rv.img, na.rm = T)/(n.col * (n.row - 1)) + sum(rh.img, na.rm = T)/(n.row * (n.col - 1)))/2
  
}
