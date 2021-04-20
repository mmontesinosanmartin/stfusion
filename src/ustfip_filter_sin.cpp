// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
int filter_sin(arma::mat x, arma::vec y){
  // correlations
  arma::vec cors(x.n_cols);
  for(int i = 0; i < x.n_cols; i++){cors.row(i) = arma::cor(x.col(i),y);}
  // spectral difference
  arma::vec difs = 1 - mean(abs(x.each_col() - y), 0).as_col();
  // similarity index
  arma::vec sind = cors/sum(cors) % difs/sum(difs);
  int mxsin = sind.index_max();
  return mxsin;
}
