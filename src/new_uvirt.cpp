// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec corr_composite(arma::mat& x,
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
    out(i) = filter_sin(xng, yng);
  }
  // return result
  return out;
}

// [[Rcpp::export]]
arma::mat local_lm(arma::mat& x,
                   arma::mat& y,
                   arma::uvec& n,
                   arma::ivec& cdims,
                   int w) {
  // initialize output
  int npx = x.n_rows;
  int ntm = x.n_cols;
  int ncl = ntm + 1;
  arma::mat out(npx, ncl);
  // inner parameters
  int nrow = cdims(0);
  int ncol = cdims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // neighboring pixels
    arma::uvec inds = get_ngbs(i, w, nrow, ncol);
    // regression data
    arma::mat yng = y.rows(inds);
    arma::mat xng = x.rows(inds);
    arma::vec kns = arma::ones(xng.n_rows);
    // coefficients
    arma::mat xinp = join_horiz(xng.col(n(i)), kns);
    arma::vec coef = arma::solve(xinp, yng);
    // saving
    arma::vec outi = arma::zeros(ncl);
    outi(n(i)) = coef[0];
    outi(ncl - 1) = coef[1];
    out.row(i) = outi.as_row();
  }
  // return result
  return out;
}

// [[Rcpp::export]]
arma::vec composite(arma::mat& x,
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

// arma::uvec get_ngbs(int i, int w, int nrow, int ncol) {
//   // current row-col
//   int row = floor(i / ncol);
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
//       out(z) = (ncol * j) + i + k;
//       z++; 
//     }
//   }
//   // return output
//   return out;
// }
// 
// int filter_sin(arma::mat x, arma::vec y){
//   // correlations
//   arma::vec cors(x.n_cols);
//   for(int i = 0; i < x.n_cols; i++){cors.row(i) = arma::cor(x.col(i),y);}
//   // spectral difference
//   arma::vec difs = 1 - abs(mean(x.each_col() - y, 0).as_col());
//   // similarity index
//   arma::vec sind = cors/sum(cors) % difs/sum(difs);
//   int mxsin = sind.index_max();
//   return mxsin;
// }
