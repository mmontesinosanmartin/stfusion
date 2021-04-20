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
arma::vec get_heterogeneity(arma::uvec& clss,
                            int w,
                            int nrow,
                            int ncol){
  // initialize
  int npx = clss.n_elem;
  arma::vec out(npx);
  // for each pixel
  for(int i = 0; i < npx; i++){
    // class of this pixel
    int clsi = clss[i];
    // pixel neighbors
    arma::uvec ngbs = get_ngbs(i, w, nrow, ncol);
    int nngb = ngbs.n_elem;
    // neighbors of the same class
    int z = 0;
    for(int j = 0; j < nngb; j++){
      int clsj = clss(ngbs[j]);
      if(clsj == clsi) z++;
    }
    // compute fraction
    out[i] = (double) z / (double) nngb;
  }
  // return result
  return out;
}
