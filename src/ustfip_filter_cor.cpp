// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
int filter_cor(arma::mat x, arma::vec y){
  // correlations
  arma::mat cors(x.n_cols, 1);
  for(int i = 0; i < x.n_cols; i++){cors.row(i) = arma::cor(x.col(i),y);}
  // maximum
  arma::vec vcors = vectorise(cors);
  int mxcor = vcors.index_max();
  return mxcor;
}
