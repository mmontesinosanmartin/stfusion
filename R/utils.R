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
  rtmp[] <- NA; stack(lapply(1:nlyrs, function(i, tmp) tmp, tmp = rtmp))
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
