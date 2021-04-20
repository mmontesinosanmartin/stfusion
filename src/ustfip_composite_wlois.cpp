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

// [[Rcpp::export]]
double wcor(arma::vec x, arma::vec y, arma::vec w){
  double mx = sum(w % x) / sum(w);
  double my = sum(w % y) / sum(w);
  double sx = sum(w % pow((x - mx), 2)) / sum(w);
  double sy = sum(w % pow((y - my), 2)) / sum(w);
  double sxy = sum(w % (x -  mx) % (y - my)) / sum(w);
  return sxy / sqrt(sx * sy);
}

// [[Rcpp::export]]
int filter_wcor(arma::mat x, arma::vec y, int cntr){
  // correlations
  arma::mat cors(x.n_cols, 1);
  for(unsigned int i = 0; i < x.n_cols; i++){
    arma::vec diff = 1. / (sqrt(pow((x.col(i) - x(cntr, i)),2)) + 1e-8);
    cors.row(i) = wcor(x.col(i), y, diff);}
  // maximum
  arma::vec vcors = vectorise(cors);
  int mxcor = vcors.index_max();
  return mxcor;
}

// [[Rcpp::export]]
arma::uvec composite_wlois(arma::mat& x,
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
    arma::uvec cntr = find(inds == i);
    // data
    arma::mat xng = x.rows(inds);
    arma::mat yng = y.rows(inds);
    // best time period
    out(i) = filter_wcor(xng, yng, cntr[0]);
  }
  // return result
  return out;
}