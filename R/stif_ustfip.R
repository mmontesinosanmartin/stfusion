#' @title Applies the USTFIP spatiotemporal image fusion model
#'
#' @description Predicts the fine image on a target date \eqn{t_k} based on
#' a series of fine images around \eqn{t_k} and a coarse image captured on
#' that date.
#' 
#' Internally applies a modified version of the Fit-FC method. Following
#' Wang and Atkinson (2018), the parameters for the Fit-FC model are set
#' \code{sngb.lr = 5}, \code{sngb.wg=15}, and \code{nsim=30} to downscale
#' from Sentinel-3 images (300x300m) to Sentinel-2 (10x10m).
#' 
#' @param f.ts a series of fine images as a list of \code{Raster*} objects
#' @param c.tk a coarse image from \eqn{t_k} as a \code{Raster*} object
#' @param c.ts a series of fine images as a list of \code{Raster*} objects (optional)
#' @param sngb.lr number of neighboring pixels for the virtual image generation
#' @param sngb.wg number of neighboring pixels for the residual compensation
#' @param nsim number of similar pixels in each neighborhood
#' @param selec character. Data selection strategy, either "optimal" (default) or "lastavl", i.e. last available .
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' @param verbose whether to notify intermediate steps (default, \code{TRUE)}
#' 
#' @return the fine image predicted at \eqn{t_k} as a \code{RasterStack}
#' 
stif_ustfip <- function(f.ts,
                        c.tk,
                        sngb.lr,
                        sngb.wg,
                        nsim,
                        c.ts = NULL,
                        selec = "optimal",
                        scale = c(0,1),
                        ncores = 1,
                        verbose = TRUE){

  # ordering by date
  f.ts <- f.ts[order(abs(as.Date(names(f.ts)) - get_dates_from_layer(c.tk)[1]))]
  # basic info.:fine
  fnm <- dim(f.ts[[1]])[1:2]; ntm <- length(f.ts)
  ftm <- raster(f.ts[[1]]); ftm[] <- NA
  # basic info.:coarse
  cnm <- dim(c.tk)[1:2]; nly <- dim(c.tk)[3]
  ctm <- raster(c.tk); ctm[] <- NA
    
  # step 1: conformation (parallel)
  c.ts <- stfusion:::.check_cseries(c.ts, f.ts, ctm, ncores)
  c.ts <- stfusion:::.img_rearrange(c.ts)
  f.ts <- stfusion:::.img_rearrange(f.ts)
  if(verbose == TRUE) message("conformation completed")
  
  # step 2: linear regression (parallel)
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
    y.mat <- as.matrix(c.tk[[i]][])
    f.mat <- as.matrix(f.ts[[i]][])
    cones <- setValues(ctm, 1)
    fones <- setValues(ftm, 1)
    
    # composite
    if(selec == "optimal") c.cmp <- composite_lois(x.mat, y.mat, cnm, sngb.lr)
    if(selec == "lastavl") c.cmp <- composite_loa(x.mat, y.mat, cnm, sngb.lr)
    f.cmp <- as.matrix(resample(setValues(ctm, c.cmp), ftm, "ngb")[])
    f.ref <- setValues(ftm,composite_genr(f.mat, f.cmp))
    
    # fitting
    coeff <- stfusion:::.gen_tmp(ctm, ntm + 1)
    coeff <- setValues(coeff, local_lm(x.mat, y.mat, c.cmp, cnm, sngb.lr))
    
    # predicting
    c.hat <- calc(coeff * stack(c.ts[[i]], cones), sum, na.rm = TRUE)
    f.hat <- calc(resample(coeff, ftm, "ngb") * stack(f.ts[[i]], fones),
                  sum, na.rm = TRUE)
    
    # edge cases
    b.msk <- setValues(ctm, apply(coeff[],1,function(x){
      ifelse(all(is.na(x)),NA,1)}))
    c.hat <- c.hat * b.msk
    f.hat <- f.hat * resample(b.msk, ftm, "ngb")
    
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
  c.virt <- .collect_parallel(tmp, 1, ctm, nly)
  f.virt <- .collect_parallel(tmp, 2, ftm, nly)
  f.reff <- .collect_parallel(tmp, 3, ftm, nly)
  if(verbose == TRUE) message("linear regression completed")
  
  # residual compensation
  delta.c <- c.tk - c.virt
  delta.f <- stfusion:::.raster_warp(delta.c, raster(f.virt), method = "cubic")
  f2.raw <- f.virt + delta.f
  
  # spatial weighting
  f2.raw <- trim(f2.raw)
  f.reff <- crop(f.reff, f2.raw)
  f2.pred <- stfusion:::.spatial_filtering_par(f.reff, f2.raw,
                                               raster(f.reff), nly, sngb.wg,
                                               nsim, spwgt, n = ncores) 
  
  # saving results
  names(f2.pred) <- names(c.tk)
  f2.pred <- clamp(f2.pred, lower = min(scale), upper = max(scale))
  
  if(verbose) out <- list(f2.lm = f2.raw, f2.pred = f2.pred)
  if(!verbose) out <- list(f2.pred = f2.pred)
  return(out)
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
  out[] <- do.call(cbind,lapply(tmp, function(x)x[[i]]))
  return(out)
}
