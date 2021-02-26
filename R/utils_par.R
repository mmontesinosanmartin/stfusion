.split_raster <- function(dims, n){
  
  rows <- dims[1]
  cols <- dims[2]
  jump <- ceiling(rows / n) 
  
  rwspl1 <- seq(1, rows, jump)
  rwspl2 <- rwspl1 + jump - 1
  clspl1 <- rep(1, n)
  clspl2 <- rep(cols, n)
  
  quadr <- cbind(rwspl1, rwspl2, clspl1, clspl2)
  quadr[,1:2] <- apply(quadr[,1:2, drop = FALSE], 2, clamp, 1, rows)
  
  return(quadr)
}

.buffer_chunk <- function(dims, chunks, w){
  
  dims <- rep(dims[1:2], each = 2)
  nmes <- colnames(chunks)
  mins <- grep("1", nmes)
  maxs <- grep("2", nmes)
  
  chunks[,mins] <- chunks[,mins] - w
  chunks[,maxs] <- chunks[,maxs] + w
  chunks <- lapply(1:4, function(i) clamp(chunks[,i], 1, dims[i]))
  chunks <- do.call(cbind, chunks)
  colnames(chunks) <- nmes
  
  return(chunks)
  
}

.spatial_filtering_par <- function(x, y, ftm, nly, w, nsim, spwgt, n = 1){
  
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
                     res <- spatial_filtering_par(x.chnk[], y.chnk[], dim(x.chnk), i.chnk, w, nsim, spwgt)
                     rm(x.chnk, y.chnk, i.chnk)
                     gc()
                     return(res)
                   }
  parallel::stopCluster(clustr)
  return(out)
}
