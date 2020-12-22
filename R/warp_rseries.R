#' Warps a series of (multi-band) images
#' 
#' @param rseries
#' @param tmpl
#' @param method
#' @param usegdal
#' 
#' @returns
#' 
warp_rseries <- function(x, tmpl, method, usegdal){
  
}

#' Warps a multi-band image
#' 
#' @param x
#' @param tmpl
#' @param method
#' @param usegdal
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
