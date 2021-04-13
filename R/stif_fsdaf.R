#' @title Applies the FSDAF spatiotemporal image fusion model
#'
#' @description Predicts the fine image at a time \eqn{t_2} based on a fine
#' image from a previous moment time \eqn{t_1} and two coarse images, from
#' \eqn{t_1} and \eqn{t_2}.
#' 
#' @param f1 a fine image from t1 as a \code{RasterStack}
#' @param c1 a coarse image from t1 as a \code{RasterStack}
#' @param c2 a coarse image from t2 as a \code{RasterStack}
#' @param w number that defines the radius of the neighborhood
#' @param nclas integer. Number of classes for the unmixing.
#' @param npure integer. Number of similar pixels in each neighborhood
#' @param scale the dynamic range of the image (default, \code{c(0,1)}
#' @param verbose whether to notify intermediate steps (default, \code{TRUE)}
#' 
#' @return the fine image predicted at t2 as a \code{RasterStack}
#' 
#' @references Zhu, X., Helmer, E. H., Gao, F., Liu, D., Chen, J., & Lefsky,
#' M. A. (2016). A flexible spatiotemporal method for fusing satellite images
#' with different resolutions. Remote Sensing of Environment, 172, 165-177.
#' 
stif_fsdaf <- function(
  c1,
  c2,
  f1,
  w = 20,
  k = 4,
  ktype = "iso",
  npure = 100,
  scale = c(0,1),
  big = FALSE,
  verbose = FALSE
){
  
  # basic information
  nbd <- nlayers(c1)
  fdm <- dim(f1)
  
  # image adjustment
  zoom <- max(ceiling(res(c1)/res(f1)))
  ctof <- resample(setValues(raster(c1), 1:ncell(c1)), f1, "ngb")
  f1 <- stfusion:::.raster_warp(f1, ctof, "cubic", usegdal = TRUE)
  
  # image classification
  if(verbose) message("1) classification")
  if(ktype == "iso") cls <- setValues(raster(f1), classify_isodata(f1[]*1e4, k) + 1)
  if(ktype == "kmeans") cls <- setValues(raster(f1), kmeans(f1[], k)$cluster)
  nclas <- length(unique(cls[]))
  if(verbose) message(paste0("number of classes: ", nclas))
  
  # neglecting extreme values in f1
  if(verbose) message("2) unmixing spectral endmembers")
  for(i in 1:nbd){
    f1.lyr <- f1[[i]]
    f1.rng <- quantile(f1[[i]], probs = c(0.0001, 0.9999))
    px.inf <- f1[[i]] <= min(f1.rng); px.sup <- f1[[i]] >= max(f1.rng)
    f1.new <- f1[[i]]; f1.new[px.inf] <- NA; f1.new[px.sup] <- NA
    f1.new[px.inf] <- min(f1.new[],na.rm = TRUE)
    f1.new[px.sup] <- max(f1.new[],na.rm = TRUE)
    f1[[i]] <- f1.new
  }
  
  # readjusting spatial extent
  f1c <- stfusion:::.raster_warp(f1, c1, "average", usegdal = TRUE)
  c1f <- stfusion:::.raster_warp(c1, f1, "near", usegdal = TRUE)
  c2f <- stfusion:::.raster_warp(c1, f1, "near", usegdal = TRUE)
  
  # get the fractures
  cls.mat <- as.matrix(cls[])
  frc <- setValues(c1, get_fractures(cls.mat, ctof[], nclas))
  het <- setValues(cls, get_heterogeneity(cls.mat, w, fdm[1],fdm[2]))
  
  # Average change rate per class
  cdt <- c2 - c1
  min.d <- cellStats(cdt, "min") - cellStats(cdt, "sd")
  max.d <- cellStats(cdt, "max") + cellStats(cdt, "sd")
  
  # temporal rate of change
  cdt.mat <- cdt[]
  frc.mat <- frc[]
  coefs <- matrix(NA, nclas, nbd)
  for(lyi in 1:nbd){
    data.lyi <- cbind(frc.mat, cdt.mat[,lyi])
    pure <- list()
    for(ci in 1:nclas){
      # select the purest
      frc.ord <- rev(order(data.lyi[,ci], decreasing = F))
      classfv <- which(data.lyi[,ci] > 0.01)
      selpure <- min(length(classfv), npure)
      data.lyci <- data.lyi[frc.ord[1:selpure],]
      # avoid (as possible) land use changes or no change
      cdt.qua <- quantile(data.lyci[,(nclas+1)],probs = c(0.1, 0.9), type = 6)
      pure[[ci]] <- data.lyci[data.lyci[,(nclas+1)]>=cdt.qua[1]&data.lyci[,(nclas+1)]<=cdt.qua[2],]
    }
    # computing rate of change
    data <- do.call(rbind,pure)
    coefs[,lyi] <- coefficients(lm(data[,(nclas+1)] ~ data[,1:nclas]+0))
  }
  df <- setValues(f1, coefs[cls[],])
  f2.tph <- clamp(f1 + df, min(scale), max(scale))
  
  # aggregate prediction to coarse resolution
  f2c <- stfusion:::.raster_warp(f2.tph, c1, "average", usegdal = TRUE)
  
  # find minimum and maximum values allowed
  min.allow <- apply(rbind(cellStats(f2.tph, "min"),cellStats(c2, "min")), 2, min)
  max.allow <- apply(rbind(cellStats(f2.tph, "max"),cellStats(c2, "max")), 2, min)
  
  # prediction using tps
  # f2.tps <- .raster_tps(f2c, stfusion:::.gen_tmp(raster(f1), 1), big = big)
  f2.tps <- stfusion:::.raster_warp(f2c, stfusion:::.gen_tmp(raster(f1), 1), "cubic", usegdal = TRUE)
  
  # computation of residuals
  predict.change <- f2c- f1c
  real.change <- c2 - c1
  change.r <- real.change - predict.change
  
  # restrict extreme values of residuals
  w.change.tps <- f2.tps - f2.tph
  dif.change <- resample(change.r, f1, "ngb")
  w.change.tps[dif.change <= 0 & w.change.tps > 0,] <- 0
  w.change.tps[dif.change >  0 & w.change.tps < 0,] <- 0
  w.change.tps <- abs(w.change.tps)
  w.uniform <- abs(dif.change)
  
  # compute residual weighs based on heterogeneity
  w.change <- (w.change.tps * het) + (w.uniform * (1 - het)) + 1e-5
  w.chnavg <- resample(stfusion:::.raster_warp(w.change, c1, "average", usegdal = TRUE), f1, "ngb")
  w.change <- w.change / w.chnavg
  
  # restrict the residuals and normalize
  w.change <- clamp(w.change, lower = -10, upper = 10)
  w.chnavg <- resample(stfusion:::.raster_warp(w.change, c1, "average", usegdal = TRUE), f1, "ngb")
  w.change <- w.change / w.chnavg
  delta <- w.change * dif.change
  
  # restrict the prediction
  f2.hat <- f2.tph + delta
  for(i in 1:nbd){
    f2i <- f2.hat[[i]]
    f2i[f2i < min.allow[i]] <- min.allow[i]
    f2i[f2i > max.allow[i]] <- max.allow[i]
    f2.hat[[i]] <- f2i
  }
  
  # spatial filtering delta
  delta <- f2.hat - f1
  f2.hat <- fsdaf_spatial_filtering(as.matrix(c1f), as.matrix(c2f), as.matrix(f1),
                                    as.matrix(delta), dim(f1), w, w, nclas,
                                    min(scale), max(scale))
  f2.hat <- setValues(f1, f2.hat)
  f2.hat <- clamp(f2.hat, min(scale), max(scale))
  return(f2.hat)
}

.raster_tps <- function(s, tmpl, big){
  # for small rasters
  if(!big) return(stack(lapply(1:nlayers(s), function(i, r){
    xy <- data.frame(xyFromCell(r[[i]], 1:ncell(r[[i]])))
    v <- getValues(r[[i]])
    xy <- xy[!is.na(v),]
    v <- v[!is.na(v)]
    tps <- Tps(xy, v)
    interpolate(tmpl, tps)
  }, r = s)))
  # for large rasters
  if(big) return(stack(lapply(1:nlayers(s), function(i, r){
    xy <- data.frame(xyFromCell(r[[i]], 1:ncell(r[[i]])))
    v <- getValues(r[[i]])
    xy <- xy[!is.na(v),]
    v <- v[!is.na(v)]
    # tps <- bigtps(xy, v)
    tps <- fastTps(xy, v, theta = 0.018)
    interpolate(tmpl, tps)
  }, r = s)))
}

