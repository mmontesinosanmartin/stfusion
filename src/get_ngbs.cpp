// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec get_ngbs(int i, int w, int nrow, int ncol) {
  // current row-col
  int row = (int) floor(i / ncol);
  int col = i - (ncol * row);
  // window limits
  int wcs = -w; if(col-w < 0) wcs = - col;
  int wce =  w; if((col + w) > (ncol - 1)) wce = ncol - col -1;
  int wrs = -w; if(row-w < 0) wrs = - row;
  int wre =  w; if(row + w > nrow - 1) wre = nrow - row - 1;
  // initialize output
  int wln = (abs(wcs) + abs(wce) + 1) * (abs(wrs) + abs(wre) + 1);
  arma::uvec out(wln);
  // sliding indices
  size_t z = 0;
  for(int j = wrs; j <= wre; j++){
    for(int k = wcs; k <= wce; k++){
      out[z] = (ncol * j) + i + k;
      z++;
    }
  }
  // return output
  return out;
}