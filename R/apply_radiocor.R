#' @title Applies the radiometric correction to multispectral images
#' 
#' @description Computes the slope and intercept of matching pairs of 
#' images from different sensors and performs the radiometric correction.
#' Applies the correction over \code{x} to resemble those in \code{y}.
#' 
#' @details Both \code{x} and \code{y} must have the same resolution and
#' number of bands. When \code{x} are fine-scale images, transformation
#' to lower scales can be done with \link{get_coarse()}. The lists must
#' be named with the capturing date of the image.
#' 
#' @param x a list of images (\code{RasterStacks}s) to be transformed
#' @param y a list of reference images (\code{RasterStack}s)
#' @param wndw integer, spatial neighbors for the local regression. Default 2.
#' @param pars whether you want to return the parameters. Default \code{FALSE}.
#' @param verbose whether to print the processing steps. Default \code{FALSE}.
#'
#' @returns a list of corrected images
#' 

apply_radiocor <- function(x, y, wndw = 2, pars = FALSE, verbose = FALSE){
  
  # matching pairs
  x.mtch <- .match_dates(names(x), names(y))
  y.mtch <- .match_dates(names(y), names(x))
  x.mtch <- x[x.mtch]
  y.mtch <- y[y.mtch]
  
  # rearrange
  x.mtch <- stfusion:::.img_rearrange(x.mtch)
  y.mtch <- stfusion:::.img_rearrange(y.mtch)
  
  # info
  tmpl <- raster(x.mtch[[1]])
  nlyr <- length(x.mtch)
  dims <- dim(x.mtch[[1]])
  
  # coefficients
  slope <- stfusion:::.gen_tmp(raster(x.mtch[[1]]), nlyr)
  inter <- stfusion:::.gen_tmp(raster(x.mtch[[1]]), nlyr)
  if(verbose) message("estimating correction parameters...")
  for(i in 1:nlyr){
    if(verbose) message("layer ", i, " is now being processed")
    x.mat <- as.matrix(x.mtch[[i]][])
    y.mat <- as.matrix(y.mtch[[i]][])
    coefs <- radio_par(x.mat, y.mat, dims, wndw)
    slope[[i]][] <- coefs[,1]
    inter[[i]][] <- coefs[,2]
  }
  
  # correction
  if(verbose) message("applying the correction...")
  x.new <- lapply(x, function(x, slope, inter){
    slope * x + inter
  }, slope = slope, inter = inter)
  names(x.new) <- names(x)
  
  # return
  out <- list(imgs = x.new)
  if(pars){
    out$slope <- slope
    out$inter <- inter
  }
  return(out)
}

.match_dates <- function(x.dte, y.dte){
  sapply(x.dte, function(xnm, ynm)
    which(ynm == xnm), ynm = x.dte)
}
