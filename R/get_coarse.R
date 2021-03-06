#' @title Obtains the coarse version of a fine-scale image
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
#' @param verbose whether to print the processing steps. Default \code{FALSE}
#' 
#' @returns a (list of) \code{RasterStack} with the coarse-scale version of \code{fimg}
#' 
get_coarse <- function(fimg,
                       tmpl,
                       slope = NULL,
                       inter = NULL,
                       blurr = TRUE,
                       radcor = FALSE,
                       verbose = FALSE){
  
  # Check inputs
  # ============
  # check images
  is.img.series <- is.list(fimg)
  if(is.img.series){
    nimg <- length(fimg)
    imnm <- names(fimg)
    nlyr <- nlayers(fimg[[1]])
    fimg <- fimg
  } else {
    nimg <- 1
    imnm <- NULL
    nlyr <- nlayers(fimg)
    imgs <- list(); imgs[[1]] <- fimg
  }
  # check templates
  is.tmp.series <- is.list(tmpl)
  if(!is.tmp.series && verbose) {
    message("we use the same template for all")
    tmpl <- list(tmpl)
  }
  # checks for rad. cor.
  rad.cor <- !is.null(slope)&!is.null(inter)&radcor == TRUE
  if(rad.cor && verbose){
    if(length(slope) != nimg | length(inter) != nimg)
      message("slope and inter must be the same length as fimg")
      stop()
  }
  
  # Processing:
  # ===========
  out <- list()
  for(i in 1:nimg){
    # Select image
    if(verbose) message(paste0("processing image ", i))
    fimgi <- fimg[[i]]
    tmpli <- tmpl[[ifelse(is.tmp.series, i, 1)]]
    # Aggregation
    fimgi[is.na(fimgi)] <- Inf
    chati <- warp_stack(fimgi, tmpli, "average", usegdal = TRUE)
    chati[is.infinite(chati)] <- NA
    # Blurring
    if(blurr) chati[] <- apply_blur(as.matrix(chati[]), dim(chati), 1)
    # Radiometric correction
    if(rad.cor) chati <- slope[[i]] * chati + inter[[i]]
    # Save
    out[[i]] <- chati
  
  }

  # return
  names(out) <- imnm
  if(!is.img.series) out <- unlist(out)
  return(out)
  
}
