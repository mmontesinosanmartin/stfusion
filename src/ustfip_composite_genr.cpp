// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec composite_genr(arma::mat& x,
                         arma::uvec& n){
  // initialize output
  int npx = x.n_rows;
  arma::vec out(npx);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // best time period
    out(i) = x(i, n(i));
  }
  // return result
  return out;
}