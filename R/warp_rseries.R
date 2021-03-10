#' Warps a multi-band image
#' 
#' @param x  a \code{RasterStack} of a multispectral image
#' @param tmpl  a \code{RasterLayer} reference template
#' @param method the warping method (e.g., "near", "average", etc.)
#' @param usegdal whether to use gdal (default=TRUE)
#' 
#' @returns
#' 
warp_stack <- function(x, tmpl, method, usegdal){
  raster::stack(lapply(1:nlayers(x),
    function(i, x, tmpl, method, usegdal){
       as(st_warp(st_as_stars(x[[i, drop = FALSE]]),
                  st_as_stars(raster(tmpl)),
                  method = method,
                  use_gdal = usegdal),
                  "Raster")
    }, x = x, tmpl = tmpl, method = method, usegdal = usegdal))
}
