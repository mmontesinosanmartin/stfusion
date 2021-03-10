// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

//' @title Calculate radiometric parameters
//' 
//' @description Calculates the slope and intercept to radiometrically correct
//' satellite images.
//' 
//' @details the function has to be applied for each band independently. 
//' Matrices x and y must have the pixels as rows and different dates as
//' columns. Both matrices must have the same length and NAs are not
//' allowed.
//' 
//' @param x a matrix with the original image (pixels x dates)
//' @param y a matrix with the reference image (pixels x dates)
//' @param dims a vector with the dimensions of the image (rows x columns)
//' @param w a number specifying the window radius (side 2w+1)
//' 
//' @returns a matrix with pixels x parameters (slope and intercept)
//' 
// [[Rcpp::export]]
arma::mat radio_par(arma::mat x, arma::mat y, arma::uvec dims, int w){
  
  // initialize
  int npx = y.n_rows;
  arma::mat out(npx, 2);

  // parameters
  int nrows = dims[0];
  int ncols = dims[1];
  
  // compute
  for(int i = 0; i < npx; i++){
    // neighbors
    arma::uvec ngbs = get_ngbs(i, w, nrows, ncols);
    // data
    arma::vec yraw = vectorise(y.rows(ngbs));
    arma::vec xraw = vectorise(x.rows(ngbs));
    // clean
    arma::uvec cmp = complete_obs(join_horiz(xraw,yraw), 1);
    if(!cmp.is_empty()){
      arma::vec xi = xraw.rows(cmp);
      arma::vec yi = yraw.rows(cmp);
      // fit
      arma::vec cv = cov(xi,yi);
      double vr = var(xi);
      double sl = cv(0)/vr;
      double it = mean(yi) - sl * mean(xi);
      // save
      out(i,0) = sl;
      out(i,1) = it;
    }
  }
  // return output
  return out;
}
