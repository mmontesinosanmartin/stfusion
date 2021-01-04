#' @title Applies the UVIRT spatiotemporal data fusion model
#'
#' @description Predicts the fine image at a time tk based on a series of fine images around
#' tk and a coarse image from that time.
#' 
#' @param f.ts  a series of fine images as a \code{RasterStack}
#' @param c2 coarse image from t2 as a \code{RasterBrick}
#' @param sngb.lr number of neighboring pixels for the virtual image generation
#' @param sngb.wg number of neighboring pixels for the residual compensation
#' @param nsim number of similar pixels in each neighborhood
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' 
#' @return the fine image predicted at tk as a \code{RasterStack}

stif_uvirt <- function(f.ts, c2, sngb.lr, sngb.wg, nsim, scale = c(0,1)){
  
  # settings
  sngb.lr <- max(1, sngb.lr)
  spwgt <- sngb.wg + 0.5
  
  # COARSE IMAGES
  # ==============
  f.ts <- .img_rearrange(f.ts)
  c.ts <- get_coarse(f.ts, raster(c2))
  
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
  
  # VIRTUAL IMAGE GENERATION
  #=========================
  # outputs
  # fine virtual image
  # coarse virtual image
  # regression coefficients
  # error image
  a.c <- list()
  b.c <- list()
  f.virt <- list()
  c.virt <- list()
  c.coefs <- .gen_tmp(ctm, ntm + 1)
  e.f <- list()
  
  # for each layer
  for(i in 1:nly){
    
    # communication is key
    message(paste("predicting layer", i, "..."))
    
    # data
    x.train <- c.ts[[i]]
    x.predi <- f.ts[[i]]
    y.train <- c2[[i]]
    
    # coarse coefficients
    x.tmat <- as.matrix(x.train[]);
    y.tmat <- as.matrix(y.train[]);
    c.coefs[] <- local_regression(x.tmat, y.tmat, cnm, sngb.lr)
    a.c[[i]] <- c.coefs[[1:ntm]]
    b.c[[i]] <- c.coefs[[ntm + 1]]
    
    # fine coefficients
    f.coefs <- resample(c.coefs, ftm, method = "ngb")
    a.f <- f.coefs[[1:ntm]]
    b.f <- f.coefs[[ntm + 1]]
    
    # virtual images
    if(ntm == 1) {
      c.virt[[i]] <- a.c[[i]] * x.train + b.c[[i]]
      f.virt[[i]] <- a.f * x.predi + b.f
    }else{
      c.virt[[i]] <- calc(a.c[[i]] * x.train, sum, na.rm = TRUE) + b.c[[i]]
      f.virt[[i]] <- calc(a.f * x.predi, sum, na.rm = TRUE) + b.f 
    }
  }
  
  # final result
  c.virt <- stack(c.virt)
  f.virt <- stack(f.virt)
  
  # WEIGHTING RESIDUALS
  # ====================
  # temporal change coarse
  delta.c <- c2 - c.virt
  delta.f <- as(.raster_warp(delta.c, f.virt, method = "cubic"), "Raster")
  f.err <- .gen_tmp(ftm, nly);
  f.err[] <- sp_pred(f.virt[], delta.f[], fnm, sngb.wg, nsim, spwgt)
  
  # Final estimate
  f2.h <- f.virt + f.err
  f2.h <- clamp(f2.h, lower = min(scale), upper = max(scale))
  names(f2.h) <- names(c2)
  return(f2.h)
}

.raster_warp <- function(r, ref, method, usegdal = TRUE){
  st_warp(st_as_stars(r), st_as_stars(ref), method = method, use_gdal = usegdal)
}

.split_raster <- function(dims, n){
  
  rows <- dims[1]
  cols <- dims[2]
  jump <- ceiling(rows / n) 
  
  rwspl1 <- seq(1, rows, jump)
  rwspl2 <- rwspl1 + jump - 1
  clspl1 <- rep(1, n)
  clspl2 <- rep(cols, n)
  
  quadr <- cbind(rwspl1, rwspl2, clspl1, clspl2)
  quadr[,1:2] <- apply(quadr[,1:2, drop = FALSE], 2, vclamp, 1, rows)
  
  return(quadr)
}

.buffer_chunk <- function(dims, chunks, w){
  
  dims <- rep(dims[1:2], each = 2)
  nmes <- colnames(chunks)
  mins <- grep("1", nmes)
  maxs <- grep("2", nmes)
  
  chunks[,mins] <- chunks[,mins] - w
  chunks[,maxs] <- chunks[,maxs] + w
  chunks <- lapply(1:4, function(i) vclamp(chunks[,i], 1, dims[i]))
  chunks <- do.call(cbind, chunks)
  colnames(chunks) <- nmes
  
  return(chunks)
  
}

.chunk_to_cells <- function(chnki, chnkwi, cpp = TRUE){
  
  nmes <- colnames(chunks)
  rows <- grep("rw", nmes)
  cols <- grep("cl", nmes)
  
  reli <- c(chnki[rows] - chnkwi[1], chnki[cols] - chnkwi[3]) + (!cpp * 1)
  rwcl <- expand.grid(reli[3]:reli[4],reli[1]:reli[2])[,2:1]
  ncol <- diff(chnkwi[cols]) + 1
  cids <- ncol * rwcl[,1] + rwcl[,2]
  
  return(cids)
}

