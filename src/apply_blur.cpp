// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec get_ngbs(int i, int w, int nrow, int ncol) {
  // current row-col
  int row = floor(i / ncol);
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
      out(z) = (ncol * j) + i + k;
      z++; 
    }
  }
  // return output
  return out;
}

// [[Rcpp::export]]
arma::vec gaussian_kernel(arma::vec x, arma::vec y, double sigma) {
  const double pi = 3.14159265358979323846;
  arma::vec dist = (pow(x, 2) + pow(y, 2));
  double cons = 1 / (2 * pow(sigma, 2));
  arma::vec ker = (cons/pi) * exp(-1 * dist/cons);
  return ker;
}

// [[Rcpp::export]]
arma::vec gaussian_filter(int i, arma::uvec& inds, int nrow, int ncol, double sigma){
  // reference row/col
  int row = floor(i / ncol);
  int col = i - (ncol * row);
  // pixel rows and columns
  arma::uvec urows = floor(inds / ncol);
  arma::uvec ucols = inds - (ncol * urows);
  arma::vec rows = arma::conv_to<arma::vec>::from(urows);
  arma::vec cols = arma::conv_to<arma::vec>::from(ucols);
  // compute the filter
  arma::vec gaus = gaussian_kernel(rows - row, cols - col, sigma);
  arma::vec filt = gaus/sum(gaus);
  return filt;
}

// [[Rcpp::export]]
arma::mat apply_blur(arma::mat& r, arma::uvec& rdims, double sigma) {
  // initialize output
  int npx = r.n_rows;
  int nbd = r.n_cols;
  arma::mat out(npx, nbd);
  // inner parameters
  int nrow = rdims(0);
  int ncol = rdims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // neighbors
    arma::uvec ngbs = get_ngbs(i, 1, nrow, ncol);
    arma::mat  rngb = r.rows(ngbs);
    // smooth
    arma::vec filt = gaussian_filter(i, ngbs, nrow, ncol, sigma);
    arma::mat rwgt = rngb.each_col() % filt;
    out.row(i) = sum(rwgt, 0);
  }
  return out;
}