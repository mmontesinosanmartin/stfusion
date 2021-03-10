#' @title Applies the Fit-FC spatiotemporal image fusion model
#'
#' @description Predicts the fine image at a time \eqn{t_2} based on a fine
#' image from a previous moment time \eqn{t_1} and two coarse images, from
#' \eqn{t_1} and \eqn{t_2}.
#' 
#' Applies the Fit-FC method as described in Wang and Atkinson (2018).
#' In the same study, the parameters of the Fit-FC model are set to
#' \code{sngb.lr = 5}, \code{sngb.wg=15}, and \code{nsim=30} to downscale
#' from Sentinel-3 images (300x300m) to Sentinel-2 (10x10m).
#' 
#' @param f1 a fine image from t1 as a \code{RasterStack}
#' @param c1 a coarse image from t1 as a \code{RasterStack}
#' @param c2 a coarse image from t2 as a \code{RasterStack}
#' @param sngb.lr number that defines the radius of the neighborhood
#' @param sngb.wg number of neighboring pixels for the residual compensation
#' @param nsim number of similar pixels in each neighborhood
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' @param verbose whether to notify intermediate steps (default, \code{TRUE)}
#' 
#' @return the fine image predicted at t2 as a \code{RasterStack}
#' 
#' @references Wang, Q., & Atkinson, P. M. (2018). 
#' Spatio-temporal fusion for daily Sentinel-2 images.
#' Remote Sensing of Environment, 204, 31-42.
#' 
stif_fitfc <- function(f1,
                       c1,
                       c2,
                       sngb.lr,
                       sngb.wg,
                       nsim,
                       scale = c(0,1),
                       verbose = TRUE) {
  # coarse info
  # n x m number of pixels
  # no. of layers
  # a template for coarse images
  cnm <- dim(c2)
  nly <- cnm[3]
  ctm <- raster(c2); ctm[] <- NA
  
  # initial settings
  w <- matrix(1,nrow=(2*sngb.lr+1),ncol=(2*sngb.lr+1))
  
  # step 1: regression fitting
  # regression parameters
  c2 <- resample(c2, c1, method = "ngb")
  regressions <- .local_trend(c1, c2, w)
  a <- regressions$slope
  b <- regressions$inter
  e <- regressions$error
  if(verbose == TRUE) message("linear regression completed")
  
  # step 2: residual compensation
  # predictions
  f2.hat <- .pixel_predicton(f1, a, b, e)
  names(f2.hat) <- names(c2)
  if(verbose == TRUE) message("residual compensation completed")
  
  # step 3: spatial filtering
  # parameters
  r.r <- sngb.wg
  r.rln <- 2 * r.r + 1
  nsim <- as.integer(min(r.rln^2, nsim))
  spw <- r.rln / 2
  
  # distances
  ngb.dis <- 1 + apply(expand.grid(-r.r:r.r,-r.r:r.r), 1, function(x)sqrt(sum(x^2)))/spw
  ngb.wgt <- 1/ ngb.dis
  
  # array conversion
  f1.mat <- extend(f1, c(r.r,r.r),value=Inf)[]
  f2.mat <- extend(f2.hat, c(r.r,r.r),value=Inf)[]
  
  # spatial filtering
  out.v <- spatial_filtering(f1.mat, f2.mat, dim(f1)[1:2], r.r, ngb.wgt, nsim)
  if(verbose == TRUE) message("spatial filtering completed")
  
  # saving results
  f2.pred <- stack(f2.hat)
  f2.pred[] <- out.v
  names(f2.pred) <- names(f2.hat)
  f2.pred <- clamp(f2.pred, lower = min(scale), upper = max(scale))
  return(f2.pred)
  
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
