// [[Rcpp::depends(RcppArmadillo)]]
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
// 
// // [[Rcpp::export]]
// int filter_cor(arma::mat x, arma::vec y){
//   // correlations
//   arma::mat cors(x.n_cols, 1);
//   for(int i = 0; i < x.n_cols; i++){cors.row(i) = arma::cor(x.col(i),y);}
//   // maximum
//   arma::vec vcors = vectorise(cors);
//   int mxcor = vcors.index_max();
//   return mxcor;
// }


// [[Rcpp::export]]
arma::uvec composite_lois(arma::mat& x,
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
    out(i) = filter_cor(xng, yng);
  }
  // return result
  return out;
}