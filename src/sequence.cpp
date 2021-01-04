// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec sequence(int n) {
  arma::uvec out = arma::linspace<arma::uvec>(0, n - 1, n);
  return out;
}
