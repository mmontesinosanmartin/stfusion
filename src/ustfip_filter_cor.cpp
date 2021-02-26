// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
int filter_cor(arma::mat x, arma::vec y){
  // correlations
  arma::mat cors(x.n_cols, 1);
  for(int i = 0; i < x.n_cols; i++){
    arma::vec xi = x.col(i);
    cors.row(i) = arma::cor(xi,y);
  }
  // maximum
  arma::vec vcors = vectorise(cors);
  int mxcor = cors.index_max();
  return mxcor;
}
