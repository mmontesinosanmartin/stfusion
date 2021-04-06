// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec get_heterogeneity(arma::uvec clss,
                            int w,
                            int nrow,
                            int ncol){
  // initialize
  int npx = clss.n_elem;
  arma::vec out(npx);
  // for each pixel
  for(int i = 0; i < npx; i++){
    // class of this pixel
    int clsi = clss(i);
    // pixel neighbors
    arma::uvec ngbs = get_ngbs(i, w, nrow, ncol);
    int nngb = size(ngbs)[0];
    // neighbors of the same class
    int z = 0;
    for(int j = 0; j < nngb; j++){
      int clsj = clss(ngbs(j));
      if(clsj == clsi) z++;
    }
    // compute fraction
    out(i) = (double) z / nngb;
  }
  // return result
  return out;
}
