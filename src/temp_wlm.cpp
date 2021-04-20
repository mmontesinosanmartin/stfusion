//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// // [[Rcpp::export]]
// arma::uvec get_ngbs(int i, int w, int nrow, int ncol) {
//   // current row-col
//   int row = (int) floor(i / ncol);
//   int col = i - (ncol * row);
//   // window limits
//   int wcs = -w; if(col-w < 0) wcs = - col;
//   int wce =  w; if((col + w) > (ncol - 1)) wce = ncol - col -1;
//   int wrs = -w; if(row-w < 0) wrs = - row;
//   int wre =  w; if(row + w > nrow - 1) wre = nrow - row - 1;
//   // initialize output
//   int wln = (abs(wcs) + abs(wce) + 1) * (abs(wrs) + abs(wre) + 1);
//   arma::uvec out(wln);
//   // sliding indices
//   size_t z = 0;
//   for(int j = wrs; j <= wre; j++){
//     for(int k = wcs; k <= wce; k++){
//       out[z] = (ncol * j) + i + k;
//       z++;
//     }
//   }
//   // return output
//   return out;
// }

//[[Rcpp::export]]
arma::vec get_wgts(arma::vec& x, int cntr){
  // sum of square differences
  arma::vec diff = x - x(cntr);
  arma::vec dfsq = 1. / (abs(diff) + 1e-8);
  // return normalized
  return sqrt(dfsq/sum(dfsq));
}

// [[Rcpp::export]]
arma::vec fit_wlm(arma::vec& x,
                  arma::mat& y,
                  arma::vec& wgt){
  // add constant
  arma::mat xmat = arma::join_rows(x, arma::ones(x.n_rows));
  // solve coefficients
  arma::vec beta = solve(xmat.each_col() % wgt, y.each_col() % wgt);
  return beta;
}

// [[Rcpp::export]]
arma::mat local_wlm(arma::mat& x,
                    arma::mat& y,
                    arma::uvec& n,
                    arma::ivec& cdims,
                    int w) {

  // initialize output
  int npx = x.n_rows;
  int ntm = x.n_cols;
  int ncl = ntm + 1;
  arma::mat out(npx, ncl);
  out.fill(arma::datum::nan);
  // inner parameters
  int nrow = cdims(0);
  int ncol = cdims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // neighboring pixels
    arma::uvec inds = get_ngbs(i, w, nrow, ncol);
    arma::uvec cntr = find(inds == i);
    // regression data
    arma::mat yng = y.rows(inds);
    arma::mat xng = x.rows(inds);
    arma::vec xin = xng.col(n(i));
    // computing weights
    arma::vec wgt = get_wgts(xin, cntr[0]);
    // when available
    if(xin.is_finite()){
      // fit
      arma::vec coef = fit_wlm(xin, yng, wgt);
      // saving
      arma::vec outi = arma::zeros(ncl);
      outi(n(i)) = coef[0];
      outi(ncl - 1) = coef[1];
      out.row(i) = outi.as_row();
    }
  }
  // return result
  return out;
}
