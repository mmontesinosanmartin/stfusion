// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat get_fractures(arma::uvec& img,
                        arma::uvec& ctof,
                        int k){
  
  arma::uvec coarse = unique(ctof);
  int ncoarse = coarse.n_elem;
  arma::mat out(ncoarse, k);
  // for each coarse pixel
  for(int i = 0; i < ncoarse; i++){
    // fine pixels within the coarse
    arma::uvec indxi = find(ctof == (i + 1));
    arma::uvec clssi = img(indxi);
    // class-based frequencies
    arma::rowvec freqs(k);
    // count frequencies
    for(int j = 0; j < k; j++){
      arma::uvec clssj = find(clssi == j + 1);
      freqs(j) = (double) size(clssj)[0];
    }
    // turn into fractions and save
    out.row(i) = freqs/indxi.n_elem ;
  }
  // return output
  return out;
}