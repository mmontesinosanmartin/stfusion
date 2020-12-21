#' Obtains the radiometric correction parameters
#' 
#' The parameters are local slopes and intercepts from regressions.
#' 
#' @param cimg coarse-scale image, as a \code{RasterStack}
#' @param chat estimated coarse-scale image, as a \code{RasterStack}
#' @param w    integer, spatial neighbors for the local regression. Default 2.
#' @param pred bool, whether you want to return the prediction. Default \code{TRUE}
#' 
#' @returns a \code{RasterStack} with the corrected image
#' 
radio_par <- function(cimg, chat, w = 2, pred = TRUE){
  
  # fit
  wmat <- matrix(1,nrow=(2*w+1),ncol=(2*w+1))
  m <- .local_trend(cimg, chat, wmat)
  a <- m$slope
  b <- m$inter
  e <- m$error; e[] <- 0
    
  # output
  if(pred) {
    out <- .pixel_predicton(chat, a, b, e)
  } else {
    out <- list(slope = stack(a), inter = stack(b), error = stack(e))
  }
  
  # return
  return(out)
  
}

#' Gets the coarse version of a fine-scale image
#' 
#' Aggregates to the coarse resolution and applies a gaussian kernel
#' as the point spread function 
#' 
#' @param fimg  fine scale image as a \code{RasterStack}
#' @param ctm   coarse image template as a \code{RasterLayer}
#' @param slope slopes image as a \code{RasterStack}
#' @param inter intercept image as a \code{RasterStack}
#' @param blur  bool, whether to apply the Gaussian kernel. Default \code{TRUE}
#' @param corr  bool, whether to apply the radiometric correction. Default \code{FALSE}
#' 
#' @returns a \code{RasterStack} with the coarse-scale version of \code{fimg}
#' 
get_coarse <- function(fimg, ctm, slope = NULL, inter = NULL, blurr = TRUE, corr = FALSE){
  
  # initialize
  nly <- nlayers(fimg)
  out <- list()
  
  # for each layer
  for(i in 1:nly){
    
    # select
    f.i <- fimg[[i, drop = FALSE]]

    # aggregation
    f.i[is.na(f.i)] <- Inf
    chat <- .raster_warp(f.i, ctm, "average")
    chat[is.infinite(chat)] <- NA
    
    # blurring
    if(blurr) chat[] <- blur(as.matrix(chat[]), dim(chat), 1)
    
    # radiometric correction
    if(!is.null(slope)&!is.null(inter)&corr == TRUE){
      a <- slope[[i, drop = FALSE]]
      b <- inter[[i, drop = FALSE]]
      e <- ctm; e[] <- 0
      chat <- .pixel_predicton(chat, a, b, e)
    }
    
    # save
    out[[i]] <- chat
  }
  
  # return
  out <- stack(out)
  names(out) <- names(fimg)
  return(out)

}

# fits local linear trends
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

# warps with stars from raster
.raster_warp <- function(r, ref, method, usegdal = TRUE){
  raster::stack(lapply(1:nlayers(r),
                       function(i, r, ref, method, usegdal){
                         as(st_warp(st_as_stars(r[[i, drop = FALSE]]),
                                    st_as_stars(raster(ref)),
                                    method = method,
                                    use_gdal = usegdal),
                            "Raster")
                       }, r = r, ref = ref, method = method, usegdal = usegdal))
}

# performs linear regression predictions
.pixel_predicton <- function(x, slope, inter, error){
  .raster_warp(slope, x, method = "near") * x +
    .raster_warp(inter, x, method = "near") +
    .raster_warp(error, x, method = "cubic")
}
