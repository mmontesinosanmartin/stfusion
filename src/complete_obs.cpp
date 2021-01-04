// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Finds rows/columns without NAs
//' 
//' @description provides the index for the rows and columns without any
//' missing value
//' 
//' @param m a dataset, as a matrix
//' @param dimn can be 0 or 1 to find NAs across rows and columns respectively
//'  
// [[Rcpp::export]]
arma::uvec complete_obs(arma::mat m, int dimn){
  
  // value replacement
  m.elem(find_finite(m)).fill(0);
  m.replace(arma::datum::nan, 1);
  // number of NAs
  arma::vec elems = sum(m, dimn).as_col();
  // not NAs
  arma::uvec comp = find(elems <= 0);
  return comp;

}