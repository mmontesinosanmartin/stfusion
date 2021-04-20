#' @title Obtains the coarse version of a fine-scale image (parallel version)
#' 
#' @description Aggregates to the coarse resolution, applies a Gaussian kernel
#' as a representation of the point spread function , and applies a
#' linear model to radiometrically correct the images.
#' 
#' @details The slope and intercept for the radiometric corrections must have the
#' same spatial resolution as the coarse template and the same number of 
#' bands as the multi-band fine scale image.
#' 
#' @param fimg  a (list of) \code{RasterStack} of the multi-band fine-scale image(s)
#' @param tmpl  a (list of) \code{RasterLayer} of template of the coarse-scale image
#' @param slope a (list of) \code{RasterStack} of slopes for the radiometric correction
#' @param inter a (list of) \code{RasterStack} of intercepts for the radiometric correction
#' @param blurr whether to apply the Gaussian kernel. Default \code{TRUE}
#' @param radcor whether to apply the radiometric correction. Default \code{FALSE}
#' @param ncores number of cores to process in parallel
#' @param verbose whether to print the processing steps. Default \code{FALSE}
#' 
#' @returns a (list of) \code{RasterStack} with the coarse-scale version of \code{fimg}
#' 
get_coarse_par <- function(fimg,
                           tmpl,
                           slope = NULL,
                           inter = NULL,
                           blurr = TRUE,
                           radcor = FALSE,
                           ncores = 1,
                           verbose = FALSE){
  
  # Check inputs
  # ============
  # check images
  is.img.series <- is.list(fimg)
  if(is.img.series){
    nimg <- length(fimg)
    imnm <- names(fimg)
    nlyr <- nlayers(fimg[[1]])
    fimgs <- fimg
  } else {
    nimg <- 1
    imnm <- NULL
    nlyr <- nlayers(fimg)
    fimgs <- list(); fimgs[[1]] <- fimg
  }
  # check templates
  is.tmp.series <- is.list(tmpl)
  if(!is.tmp.series && verbose) {
    message("we use the same template for all")
    tmpl <- list(tmpl)
  }
  # checks for rad. cor.
  rad.cor <- !is.null(slope)&!is.null(inter)
  if(rad.cor && verbose){
    if(length(slope) != nimg | length(inter) != nimg)
      message("slope and inter must be the same length as fimg")
      stop()
  }
  
  # Processing:
  # ===========
  clustr <- makeCluster(ncores)
  doParallel::registerDoParallel(clustr)
  out <- foreach(i = 1:nimg, .packages = c("raster", "stfusion")) %dopar% {
    # Select image
    if(verbose) message(paste0("processing image ", i))
    fimgi <- fimgs[[i]]
    # fimgi <- disaggregate(fimg[[i]], 1, method = "")
    tmpli <- tmpl[[ifelse(is.tmp.series, i, 1)]]
    # Aggregation
    fimgi[is.na(fimgi)] <- Inf
    chati <- warp_stack(fimgi, tmpli, "average", usegdal = TRUE)
    chati[is.infinite(chati)] <- NA
    # Blurring
    if(blurr) chati[] <- apply_blur(as.matrix(chati[]), dim(chati), 1)
    # Radiometric correction
    if(rad.cor) chati <- slope[[i]] * chati + inter[[i]]
    # return
    (chati)
  }
  parallel::stopCluster(clustr)
  
  # return
  names(out) <- imnm
  if(!is.img.series) out <- unlist(out)
  return(out)
  
}
