// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
int filter_loa(arma::mat x){
  // finite columns?
  int ntk = x.n_cols;
  int loa = ntk;
  // closest layer with finite values
  // assumes that when i grows the
  // distance from current day also
  // grows
  for(int i = 0; i < ntk; i++){
    if(x.col(i).is_finite()){
      loa = i;
      break;
    }
  }
  // return
  return loa;
}

// [[Rcpp::export]]
arma::uvec composite_loa(arma::mat& x,
                         arma::vec& y,
                         arma::ivec& cdims,
                         int w){
  // initialize output
  int npx = x.n_rows;
  arma::uvec out(npx);
  // inner parameters
  int nrow = cdims(0);
  int ncol = cdims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // neighbors
    arma::uvec inds = get_ngbs(i, w, nrow, ncol);
    // data
    arma::mat xng = x.rows(inds);
    arma::mat yng = y.rows(inds);
    // best time period
    out(i) = filter_loa(xng);
  }
  // return result
  return out;
}