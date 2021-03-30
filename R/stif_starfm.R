#' @title Applies the STARFM image fusion model
#' 
#' @description Predicts the fine image at a time \eqn{t_2} based on a fine
#' image from a previous moment time \eqn{t_1} and two coarse images, from
#' \eqn{t_1} and \eqn{t_2}.
#' 
#' This is a fast implementation of the STARFM algorithm for the specific
#' situation in which a single pair of images is used to predict the fine
#' scale image at \eqn{t_2}.
#' 
#' @param c1 a coarse image from t1 as a \code{RasterStack}
#' @param c2 a coarse image from t2 as a \code{RasterStack}
#' @param f1 a fine image from t1 as a \code{RasterStack}
#' @param red.band integer. The number of the red band (default, 1)
#' @param nir.band integer. The number of the NIR band (default, 4)
#' @param w integer. The radius of the neighboring area (default, 20)
#' @param ns integer. Number of similar pixels considered for the spatial weighting (default, 20)
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' @param verbose whether to notify intermediate steps (default, \code{TRUE)}
#' 
#' @return the fine image predicted at t2 as a \code{RasterStack}
#' 
#' @references Gao, F., Masek, J., Schwaller, M., & Hall, F. (2006).
#' On the blending of the Landsat and MODIS surface reflectance: Predicting
#' daily Landsat surface reflectance. IEEE Transactions on Geoscience and
#' Remote sensing, 44(8), 2207-2218.
#'  
stif_starfm <- function(c1,
                        c2,
                        f1,
                        red.band = 1,
                        nir.band = 4,
                        w  = 20,
                        ns = 20,
                        scale = c(0,1),
                        verbose = TRUE){
  
  # re-sampling images
  no.bnds <- dim(c1)[3]
  c1.intp <- .raster_warp(c1, f1, "cubic", use_gdal = TRUE)
  c2.intp <- .raster_warp(c2, f1, "cubic", use_gdal = TRUE)
  
  # spatial weighting
  out <- .gen_tmpl(raster(f1), nlayers(f1))
  out[] <- starfm_spatial_filter(c1 = c1.intp[], c2 = c2.intp[], f1 = f1[],
                                 dims = dim(f1), redb = red.band,
                                 nirb = nir.band, w, ns)
  
  # changing format
  out <- clamp(out, min(scale), max(scale))
  names(out) <- names(c2)
  return(out)
}

.raster_warp <- function(r, ref, method, use_gdal){
  as(st_warp(st_as_stars(r), st_as_stars(ref), method, use_gdal = use_gdal), "Raster")
}

.gen_tmpl <- function(tmpl, nlyr){
  tmpl[] <- NA
  stack(lapply(1:nlyr, function(i, r) r, r = tmpl))
}