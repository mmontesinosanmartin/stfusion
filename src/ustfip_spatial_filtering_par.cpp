// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat spatial_filtering_par(arma::mat ref,
                                arma::mat img,
                                arma::uvec dims,
                                arma::uvec indx,
                                int w,
                                int nng,
                                double spwgt) {
  
  // initialize
  int npx = size(indx)[0];
  int nbd = ref.n_cols;
  arma::mat out(npx, nbd);

  // parameters
  int nrows = dims(0);
  int ncols = dims(1);
  arma::uvec bst = arma::linspace<arma::uvec>(0, nng - 1, nng);
  
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // stop if needed
    Rcpp::checkUserInterrupt();
    // current row-col
    int cll = indx[i];
    int row = floor(cll / ncols);
    int col = cll - (ncols * row);
    // window limits
    int wcs = -w; if((col-w) < 0) wcs = - col;
    int wce =  w; if((col+w) > (ncols - 1)) wce = ncols - col - 1;
    int wrs = -w; if((row-w) < 0) wrs = - row;
    int wre =  w; if((row+w) > (nrows - 1)) wre = nrows - row - 1;
    // initialize output
    int wln = (abs(wcs) + abs(wce) + 1) * (abs(wrs) + abs(wre) + 1);
    arma::vec simis(wln);
    arma::vec indst(wln);
    arma::mat imng(wln,nbd);
    // sliding indexes
    size_t z = 0;
    for(int j = wrs; j <= wre; j++){
      for(int k = wcs; k <= wce; k++){
        // similarity, distance, and neighbors
        int ind = (ncols * j) + cll + k;
        simis(z) = sum(pow((ref.row(ind) - ref.row(cll)), 2));
        indst(z) = 1. / (1 + sqrt(pow(j, 2) + pow(k, 2)) / spwgt);
        imng.row(z) = img.row(ind);
        // next
        z++; 
      }
    }
    
    // most similar
    arma::uvec sind = arma::sort_index(simis,"ascending");
    arma::uvec bstp = sind(bst);
    arma::mat  ibst = imng.rows(bstp);
    // normalize weights
    arma::vec  wgts = indst(bstp) / (sum(indst(bstp)) + 1e-8);
    // prediction
    arma::rowvec  pred = sum(ibst.each_col() % wgts, 0);
    out.row(i) = pred;
  }
  //result
  return out;
}