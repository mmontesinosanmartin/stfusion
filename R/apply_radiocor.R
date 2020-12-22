#' Applies the radiometric correction to a series of images
#' 
#' Calculates and applies the radiometric corrections to a series of
#' images captured on the same dates.
#' 
#' The parameters are local slopes and intercepts from regressions.
#' 
#' @param cimg coarse-scale images, as a \code{RasterStack}
#' @param chat estimated coarse-scale image, as a \code{RasterStack}
#' @param wndw integer, spatial neighbors for the local regression. Default 2.
#' @param pars whether you want to return the parameters. Default \code{FALSE}.
#' @param verbose whether to print the processing steps. Default \code{FALSE}.
#' 
#' @returns a \code{RasterStack} with the corrected image
#' 
apply_radiocor <- function(cimg, chat, wndw = 2, pars = FALSE, verbose = FALSE){
  
  if(!is.list(cimg)){
    message("cimg must be a list")
    stop()
  }
  if(!is.list(chat)){
    message("chat must be a list")
    stop()
  }
  nimg <- length(cimg)
  if(nimg != length(chat)){
    message("cimg and chat must be the same length")
    stop()
  }
  
  out <- list()
  for(i in 1:nimg){
    if(verbose) message(paste0("processing image ", i))
    out[[i]] <- radiocor(cimg[[i]], chat[[i]], wndw = wndw, pars = pars)
  }
  return(out)
}

radiocor <- function(cimg, chat, wndw = wndw, pars = FALSE){
  
  wmat <- matrix(1,nrow=(2*wndw+1),ncol=(2*wndw+1))
  m <- .local_trend(cimg, chat, wmat)
  if(pars) {
    out <- list(slope = stack(a), inter = stack(b), error = stack(e))
  } else {
    out <- m$slope * chat + m$inter
  }
  return(out)
  
}

.local_trend <- function(c1, c2, w){
  
  c1.avg <- stack(lapply(as.list(c1),function(x, w, fun, na.rm, pad, padValue){
    focal(x = x, w = w, fun = fun, na.rm = na.rm, pad = pad, padValue = padValue)
  }, w = w, fun = mean, na.rm = TRUE, pad = TRUE, padValue = NA))
  c2.avg <- stack(lapply(as.list(c2),function(x, w, fun, na.rm, pad, padValue){
    focal(x = x, w = w, fun = fun, na.rm = na.rm, pad = pad, padValue = padValue)
  }, w = w, fun = mean, na.rm = TRUE, pad = TRUE, padValue = NA))
  cx.cov <- stack(lapply(as.list(c1 * c2),function(x, w, fun, na.rm, pad, padValue){
    focal(x = x, w = w, fun = fun, na.rm = na.rm, pad = pad, padValue = padValue)
  }, w = w, fun = mean, na.rm = TRUE, pad = TRUE, padValue = NA)) - (c1.avg * c2.avg)
  c1.var <- stack(lapply(as.list(c1 * c1),function(x, w, fun, na.rm, pad, padValue){
    focal(x = x, w = w, fun = fun, na.rm = na.rm, pad = pad, padValue = padValue)
  }, w = w, fun = mean, na.rm = TRUE, pad = TRUE, padValue = NA)) - (c1.avg * c1.avg)
  
  a <- cx.cov / c1.var
  b <- c2.avg - a * c1.avg
  e <- c2 - a * c1 - b
  
  return(list(
    slope = a,
    inter = b,
    error = e
  ))
}