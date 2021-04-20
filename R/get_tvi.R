#' Calculates the Normalized Temporal Variation Index (NTVI)
#' 
#' Quantifies the temporal variation of a pair of images
#' 
#' Both \code{s1} and \code{s2} are assumed to be images of the same size.
#' 
#' The TVI ranges between 0 and 1, where 0 means no change and 1 represents
#' a large temporal variation. Temporal variations are evaluated locally in
#' spatial neighborhoods (see Thonfeld et al, 2016).  Considering a
#' neighborhood avoids geometric distortions and missregistration issues.The
#' size of the neighborhood is (2w + 1) on the horizontal and vertical axes.
#' Default value is set to 16, which is the ratio between Landsat and MODIS
#' pixels.
#' 
#' @param s1 a \code{RasterStack} at the time t1. See details.
#' @param s2 a \code{RasterStack} at the time t2 (where t2 > t1). See details.
#' @param w a number specifying the size of an spatial neighborhood. Default to
#' 16 pixels. See details.
#' @param scale the dynamic range of the image (default, \code{c(0,255)})
#' 
#' @returns \code{RasterStack} with the temporal variation at each pixel
#' 
#' @references Thonfeld, F., Feilhauer, H., Braun, M., & Menz, G. (2016).
#' Robust Change Vector Analysis (RCVA) for multi-sensor very high resolution
#' optical satellite data. International journal of applied earth observation
#' and geoinformation, 50, 131-140.
#' 
#' @example
#' library(raster)
#' r1 <- r2 <- raster(nrow = 100, ncol = 100)
#' r1[] <- rnorm(100^2,100)
#' r2[] <- rnorm(100^2,100)
#' s1 <- stack(r1,r1)
#' s2 <- stack(r2,r2)
#' get_tvi(s1,s2,w=5) 
#' 
# get_tvi <- function(s1, s2, w = 1){
# 
#   # input re-format
#   w <- as.integer(w)
#   sdim <- dim(f1)[1:2]
#   m1 <- extend(s1, c(w,w), value=Inf)[]
#   m2 <- extend(s2, c(w,w), value=Inf)[]
#   
#   # output
#   out <- raster(s1)
#   
#   # apply algorithm
#   out[] <- tvar_cpp(m1, m2, sdim, w)
#   
#   return(out)
# }
# 
# cppFunction(depends="RcppArmadillo",
# 'arma::vec tvar_cpp(arma::mat s1,
#                     arma::mat s2,
#                     arma::ivec sdim,
#                     int w) {
# 
#  int clp = sdim[1] + (2 * w);
#  int cmx = sdim[1] + w;
#  int rwp = sdim[0] + (2 * w);
#  int str = clp * w;
#  int end = clp * (rwp - w);
# 
#  arma::vec diffa(s1.n_cols);
#  arma::vec diffb(s1.n_cols);
#  arma::vec out(prod(sdim));
#  
#  // for each pixel
#  for(int i = str; i < end; i++){
# 
#   // select the neighborhood
#   int col = i  % clp;
#   int row = floor(i /clp);
#   
#   if((col >= w) & (col < cmx)) {
#    
#    diffa.fill(999);
#    diffb.fill(999);
#    
#    for(int j = -w; j <= w; j++) {
#     for(int k = -w; k <= w; k++) {
#      int ind = (clp * j) + i + k;
#      arma::rowvec diffai = clamp(s2.row(ind) - s1.row(i),0,1);
#      arma::rowvec diffbi = clamp(s1.row(ind) - s2.row(i),0,1);
#      for(int l = 0; l < s2.n_cols; l++){
#        if(diffai(l) < diffa(l)) {
#          diffa(l) = diffai(l);
#        }
#        if(diffbi(l) < diffb(l)) {
#          diffb(l) = diffbi(l);
#        }
#      }
#     }
#    }
#    
#    arma::uvec diffs0 = find(diffa == 0);
#    if(diffs0.n_elem > 0) {
#      diffa(diffs0) = diffb(diffs0) * -1;
#    }
#    
#    int outi = (row-w) * sdim[1] + (col - w);
#    arma::vec diffsq = square(diffa);
#    out(outi) = sqrt(sum(diffsq));
# 
#   }
#  }
# 
#  return(out);
# }')
