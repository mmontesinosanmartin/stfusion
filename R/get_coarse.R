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
#' @param blur  whether to apply the Gaussian kernel. Default \code{TRUE}
#' @param corr  whether to apply the radiometric correction. Default \code{FALSE}
#' @param verbose whether to print the processing steps. Default \code{FALSE}
#' 
#' @returns a (list of) \code{RasterStack} with the coarse-scale version of \code{fimg}
#' 
get_coarse <- function(fimg,
                       tmpl,
                       slope = NULL,
                       inter = NULL,
                       blurr = TRUE,
                       corr = FALSE,
                       verbose = FALSE){
  
  # check imgages
  is.img.series <- is.list(fimg)
  if(is.img.series){
    nimg <- length(fimg)
    nlyr <- nlayers(fimg[[1]])
    fimg <- fimg
  } else {
    nimg <- 1
    nlyr <- nlayers(fimg)
    imgs <- list(); imgs[[1]] <- fimg
  }
  
  # check templates
  is.tmp.series <- is.list(tmpl)
  if(!is.tmp.series) {
    message("we use the same template for all")
    tmpl <- list(tmpl)
  }

  # checks for rad. cor.
  rad.cor <- !is.null(slope)&!is.null(inter)&corr == TRUE
  if(rad.cor){
    if(length(slope) != nimg | length(inter) != nimg)
      message("slope and inter must be the same length as fimg")
      stop()
  }
  
  # Processing:
  # Select image and layer
  # Aggregate to the coarse-scale resolution (ignore clouds)
  # Blur using a Gaussian kernel (sigma = 1)
  # Apply the radiometric correction
  out <- list()
  for(i in 1:nimg){
    
    if(verbose)message(paste0("processing image ", i))
    
    fimgi <- fimg[[i]]
    tmpli <- tmpl[[ifelse(is.tmp.series, i, 1)]]
    
    fimgi[is.na(fimgi)] <- Inf
    chati <- warp_stack(fimgi, tmpli, "average", usegdal = TRUE)
    chati[is.infinite(chati)] <- NA
    
    if(blurr) chati[] <- apply_blur(as.matrix(chati[]), dim(chati), 1)
    if(rad.cor) chati <- slope[[i]] * chati + inter[[i]]
    
    out[[i]] <- chati
  
  }

  # return
  if(!is.img.series) out <- unlist(out)
  return(out)
  
}
