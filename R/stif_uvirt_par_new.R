stif_uvirt_par_new <- function(f.ts,
                               c.ts = NULL,
                               c2,
                               sngb.lr,
                               sngb.wg,
                               nsim,
                               scale = c(0,1),
                               ncores = 1){

  # BASIC INFO
  # ==========
  # fine
  fnm <- dim(f.ts[[1]])[1:2]; ntm <- length(f.ts)
  ftm <- raster(f.ts[[1]]); ftm[] <- NA
  # coarse
  cnm <- dim(c2)[1:2]; nly <- dim(c2)[3]
  ctm <- raster(c2); ctm[] <- NA
    
  # COARSE IMAGE GENERATION
  # =======================
  c.ts <- .check_cseries(c.ts, f.ts, ctm, ncores)
  c.ts <- stfusion:::.img_rearrange(c.ts)
  f.ts <- stfusion:::.img_rearrange(f.ts)
  
  # LINEAR REGRESSION
  # =================
  # parameters
  sngb.lr <- max(1, sngb.lr)
  spwgt <- sngb.wg + 0.5
  # open
  clustr <- makeCluster(ncores)
  doParallel::registerDoParallel(clustr)
  
  # layer-wise
  tmp <- foreach(i = 1:nly, .packages = c("raster", "stfusion")) %dopar% {
    
    # basic info
    x.mat <- as.matrix(c.ts[[i]][])
    y.mat <- as.matrix(c2[[i]][])
    f.mat <- as.matrix(f.ts[[i]][])
    cones <- setValues(ctm, 1)
    fones <- setValues(ftm, 1)
    
    # composite
    c.cmp <- corr_composite(x.mat, y.mat, cnm, sngb.lr)
    f.cmp <- as.matrix(resample(setValues(ctm, c.cmp), ftm, "ngb")[])
    f.ref <- setValues(ftm,composite(f.mat, f.cmp))
    
    # fitting
    coeff <- stfusion:::.gen_tmp(ctm, ntm + 1)
    coeff <- setValues(coeff, local_lm(x.mat,y.mat, c.cmp, cnm, sngb.lr))
    
    # predicting
    c.hat <- calc(coeff * stack(c.ts[[i]], cones), sum, na.rm = TRUE)
    f.hat <- calc(resample(coeff, ftm, "ngb") * stack(f.ts[[i]], fones), sum, na.rm = TRUE)
    
    # saving
    out <- list(c.hat[], f.hat[], f.ref[])
    rm(x.mat, y.mat, fones, cones,
       c.cmp, f.cmp, f.ref,
       coeff, c.hat, f.hat)
    gc()
    return(out)
  }
  
  # close
  parallel::stopCluster(clustr) 
  
  # collect
  c.virt <- stfusion:::.collect_parallel(tmp, 1, ctm, nly)
  f.virt <- stfusion:::.collect_parallel(tmp, 1, ctm, nly)
  f.reff <- stfusion:::.collect_parallel(tmp, 1, ctm, nly)
  
  # WEIGHTING RESIDUALS
  # ====================
  delta.c <- c2 - c.virt
  delta.f <- stfusion:::.raster_warp(delta.c, raster(f.virt), method = "cubic")
  f2.raw <- f.virt + delta.f
  f2.hat <- stfusion:::.sp_pred_par(f.reff, f2.raw, ftm, nly, sngb.wg, nsim, spwgt, n = ncores)  
  f2.hat <- clamp(f2.hat, min(scale), max(scale))
  return(f2.hat)
}

# Helper: check/generate coarse images when needed
.check_cseries <- function(c.ts, f.ts, ctm, ncores){
  # if there is some coarse
  if(!is.null(c.ts)){
    # dates with missing coarse
    f.dates <- names(f.ts)
    c.dates <- names(c.ts)
    miss <- which(!(f.dates %in% c.dates))
    # combine existing with generated
    if(length(miss) > 0){
      c.ts <- c(c.ts, get_coarse_par(f.ts[miss], ctm, ncores = ncores))
    }
  # if there is none
  } else {
    c.ts <- get_coarse_par(f.ts, ctm, ncores = ncores)
  }
  # ensure same order
  c.ts <- c.ts[names(f.ts)]
  return(c.ts)
}

# Helper: collect results from parallel execution
.collect_parallel <- function(tmp, i, tmpl, nly){
  out <- stfusion:::.gen_tmp(tmpl, nly)
  out <- do.call(cbind,lapply(tmp, function(x)x[[i]]))
  return(out)
}
