#' @title Applies the UVIRT spatio-temporal data fusion model
#'
#' @description Predicts the fine image at a time tk based on a series of fine images around
#' tk and a coarse image from that time.
#' 
#' @param f.ts a series of fine images as a list of \code{RasterStack}s
#' @param c2 coarse image at tk as a \code{RasterStack}
#' @param sngb.lr number of neighboring pixels for the virtual image generation
#' @param sngb.wg number of neighboring pixels for the residual compensation
#' @param nsim number of similar pixels in each neighborhood
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' @param ncores number of cores for parallel processing
#' 
#' @return the fine image predicted at tk as a \code{RasterStack}

stif_uvirt_par <- function(f.ts, c2, sngb.lr, sngb.wg, nsim, scale = c(0,1), ncores=1){
  
  # settings
  sngb.lr <- max(1, sngb.lr)
  spwgt <- sngb.wg + 0.5
  
  # COARSE IMAGES
  # ==============
  f.ts <- stfusion:::.img_rearrange(f.ts)
  c.ts <- get_coarse_par(f.ts, raster(c2), ncores = ncores)
  
  # VIRTUAL IMAGE GENERATION
  # ========================
  # fine info
  # n x m number of pixels
  # no. images in the time series
  # a template for fine images
  fnm <- dim(f.ts[[1]])
  ntm <- fnm[3]
  ftm <- raster(f.ts[[1]]); ftm[] <- NA
  
  # coarse info
  # n x m number of pixels
  # no. of layers
  # a template for coarse images
  cnm <- dim(c2)
  nly <- cnm[3]
  ctm <- raster(c2); ctm[] <- NA
  
  # outputs
  # fine virtual image
  # coarse virtual image
  # regression coefficients
  # error image
  c.coefs <- stfusion:::.gen_tmp(ctm, ntm + 1)
  clustr <- makeCluster(ncores)
  doParallel::registerDoParallel(clustr)

  # for each layer
  tmp <- foreach(i = 1:nly, .packages = c("raster", "stfusion")) %dopar% {

    # data
    x.train <- c.ts[[i]]
    x.predi <- f.ts[[i]]
    y.train <- c2[[i]]
    
    # coarse coefficients
    x.tmat <- as.matrix(x.train[]);
    y.tmat <- as.matrix(y.train[]);
    c.coefs[] <- local_regression(x.tmat, y.tmat, cnm, sngb.lr)
    a.c <- c.coefs[[1:ntm]]
    b.c <- c.coefs[[ntm + 1]]
    
    # fine coefficients
    f.coefs <- resample(c.coefs, ftm, method = "ngb")
    a.f <- f.coefs[[1:ntm]]
    b.f <- f.coefs[[ntm + 1]]
    
    # virtual images
    if(ntm == 1) {
      c.virt <- a.c * x.train + b.c
      f.virt <- a.f * x.predi + b.f
    }else{
      c.virt <- calc(a.c * x.train, sum, na.rm = TRUE) + b.c
      f.virt <- calc(a.f * x.predi, sum, na.rm = TRUE) + b.f 
    }
    
    # housekeeping
    out <- list(c.virt, f.virt)
    rm(x.train, x.predi, y.train,
       x.tmat, y.tmat, a.c, b.c,
       f.coefs, a.f, b.f, c.virt,
       f.virt)
    gc()
    return(out)
  }
  parallel::stopCluster(clustr)
  
  # final result
  c.virt <- stack(lapply(tmp, function(x)x[[1]]))
  f.virt <- stack(lapply(tmp, function(x)x[[2]]))
  
  # WEIGHTING RESIDUALS
  # ====================
  # temporal change coarse
  delta.c <- c2 - c.virt
  delta.f <- .raster_warp(delta.c, raster(f.virt), method = "cubic")
  f.err <- .sp_pred_par(f.virt, delta.f,
                        ftm, nly,
                        sngb.wg, nsim, spwgt,
                        n = ncores)
  
  # Final estimate
  f2.h <- f.virt + f.err
  f2.h <- clamp(f2.h, lower = min(scale), upper = max(scale))
  names(f2.h) <- names(c2)
  return(f2.h)
}