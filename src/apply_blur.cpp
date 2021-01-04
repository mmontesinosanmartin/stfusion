// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

arma::vec gaussian_kernel(arma::vec x, arma::vec y, double sigma) {
  const double pi = 3.14159265358979323846;
  arma::vec dist = (pow(x, 2) + pow(y, 2));
  double cons = 1 / (2 * pow(sigma, 2));
  arma::vec ker = (cons/pi) * exp(-1 * dist/cons);
  return ker;
}

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

//' @title Blurs an image
//' 
//' @description blurs an image using a Gaussian kernel, which is a usual
//' way to describe the point spread function
//' 
//' @param r a matrix with pixel by bands dimensions
//' @param rdims a vector specifying the dimensions of the image
//' @param sigma standard deviation of the Gaussian kernel, usually equal to 1
//' 
//' @retruns a blurred image as a matrix
//' 
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

