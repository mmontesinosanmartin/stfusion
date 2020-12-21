#' Applies the UVIRT spatiotemporal data fusion model
#'
#' Predicts the fine image at a time tk based on a series of fine images around
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

uvirt_par <- function(f.ts, c2, sngb.lr, sngb.wg, nsim, scale = c(0,1), ncores=1){
  
  # settings
  sngb.lr <- max(1, sngb.lr)
  spwgt <- sngb.wg + 0.5
  
  # COARSE IMAGES
  # ==============
  f.ts <- .img_rearrange(f.ts)
  c.ts <- lapply(f.ts, function(x, ctm){
    x[is.na(x)] <- Inf
    c.hat <- .raster_warp(x, ctm, "average")
    c.hat[is.infinite(c.hat)] <- NA
    c.hat[] <- blur(as.matrix(c.hat[]),dim(c.hat), 1)
    c.hat
  },ctm = raster(c2))
  
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
  c.coefs <- .gen_tmp(ctm, ntm + 1)
  clustr <- makeCluster(ncores)
  doParallel::registerDoParallel(clustr)

  # for each layer
  tmp <- foreach(i = 1:nly,
                   .packages = c("raster", "stfusion")) %dopar% {

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

.img_rearrange <- function(img.ls) {
  n.imgs <- length(img.ls)
  n.lyrs <- nlayers(img.ls[[1]])
  out <- list()
  for(i in 1:n.lyrs)
  {
    outi <- list()
    for(j in 1:n.imgs)
    {
      outi[[j]] <- img.ls[[j]][[i]]
    }
    out[[i]] <- stack(outi)
  }
  return(out)
}

.gen_tmp <- function(rtmp, nlyrs) {
  rtmp[] <- NA
  stack(lapply(1:nlyrs, function(i, tmp) tmp, tmp = rtmp))
}

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

.sp_pred_par <- function(x, y, ftm, nly, w, nsim, spwgt, n = 1){
  
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
  
  fnm <- dim(ftm)
  chunks <- .split_raster(fnm, n)
  chunkw <- .buffer_chunk(fnm, chunks, w)
  
  clustr <- makeCluster(n)
  doParallel::registerDoParallel(clustr)
  out <- .gen_tmp(ftm, nly)
  
  out[] <- foreach(i = 1:n,
                   .combine = 'rbind',
                   .packages = c("raster", "stfusion")) %dopar% {
                     
      x.chnk <- crop(x, extent(x,chunkw[i,1],chunkw[i,2],chunkw[i,3],chunkw[i,4]))
      y.chnk <- crop(y, extent(y,chunkw[i,1],chunkw[i,2],chunkw[i,3],chunkw[i,4]))
      i.chnk <- .chunk_to_cells(chunks[i,], chunkw[i,])
      res <- sp_pred_par(x.chnk[], y.chnk[], dim(x.chnk), i.chnk, w, nsim, spwgt)
      rm(x.chnk, y.chnk, i.chnk)
      gc()
      return(res)
  }
  parallel::stopCluster(clustr)
  return(out)
}
