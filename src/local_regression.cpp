// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

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

// [[Rcpp::export]]
int filter_sin(arma::mat x, arma::vec y){
  // correlations
  arma::vec cors(x.n_cols);
  for(int i = 0; i < x.n_cols; i++){cors.row(i) = arma::cor(x.col(i),y);}
  // spectral difference
  arma::vec difs = 1 - abs(mean(x.each_col() - y, 0).as_col());
  // similarity index
  arma::vec sind = cors/sum(cors) % difs/sum(difs);
  int mxsin = sind.index_max();
  return mxsin;
}

// [[Rcpp::export]]
arma::mat local_regression(arma::mat& x, arma::mat& y, arma::ivec& cdims, int w) {
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
  // int i = 14599;
    // neighboring pixels
    arma::uvec inds = get_ngbs(i, w, nrow, ncol);
    // regression data
    arma::mat xng = x.rows(inds);
    arma::vec kns = arma::ones(xng.n_rows);
    arma::mat yng = y.rows(inds);
    // best time period
    int mxcor = filter_sin(xng, yng);
    // coefficients
    arma::mat xinp = join_horiz(xng.col(mxcor), kns);
    arma::vec coef = arma::solve(xinp, yng);
    // saving
    arma::vec outi = arma::zeros(ncl);
    outi(mxcor) = coef[0];
    outi(ncl - 1) = coef[1];
    out.row(i) = outi.as_row();
  }
  // return result
  return out;
}