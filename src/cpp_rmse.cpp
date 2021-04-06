// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec cpp_rmse(arma::mat& x, arma::mat& y, bool byband){
  arma::rowvec rmse = sqrt(mean(pow(x-y,2),0));
  arma::rowvec out;
  if(!byband){
    out  = mean(rmse);
  } else {
    out = rmse;
  }
  return out;
}