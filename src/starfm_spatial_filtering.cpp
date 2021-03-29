// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat starfm_spatial_filter(arma::mat& c1,
                                arma::mat& c2,
                                arma::mat& f1,
                                arma::uvec dims,
                                unsigned int redb,
                                unsigned int nirb,
                                int w, int ns) {
  
  // initialize
  int npx = f1.n_rows;
  int nbd = f1.n_cols;
  arma::mat out(npx, nbd);
  
  // parameters
  int nrows = dims(0);
  int ncols = dims(1);
  int ngl = (2 * w) + 1;
  arma::uvec bst = arma::regspace<arma::uvec>(0, (ns - 1));
  arma::uvec rfb = {redb, nirb};
  
  // for each pixel
  for(int i = 0; i < npx; i++) {
  // int i = 979301;
  
    // stop if needed
    Rcpp::checkUserInterrupt();
    
    // current row-col
    int row = floor(i / ncols);
    int col = i - (ncols * row);
    // Rcout << row << " " << col << std::endl;
    // window limits
    int wcs = -w; if((col-w) < 0) wcs = - col;
    int wce =  w; if((col+w) > (ncols - 1)) wce = ncols - col - 1;
    int wrs = -w; if((row-w) < 0) wrs = - row;
    int wre =  w; if((row+w) > (nrows - 1)) wre = nrows - row - 1;
    // initialize output
    int wln = (abs(wcs) + abs(wce) + 1) * (abs(wrs) + abs(wre) + 1);
    arma::uvec inds(wln);
    arma::vec simis(wln);
    arma::vec indst(wln);
    arma::mat difng(wln,nbd);
    // sliding indexes
    size_t z = 0;
    for(int j = wcs; j <= wce; j++){
      for(int k = wrs; k <= wre; k++){
        // similarity, distance, and neighbors
        int ind = (ncols * k) + i + j;
        simis(z) = pow(f1(ind, (redb - 1)) - f1(i, (redb - 1)), 2) + pow(f1(ind, (nirb - 1)) - f1(i, (nirb - 1)), 2);
        double ds = sqrt(pow(k, 2) + pow(j, 2));
        indst(z) = 1. / (1 + (ds/(ngl * 0.5)));
        difng.row(z) = c2.row(ind) - c1.row(ind);
        // next
        z++; 
      }
    }

     // most similar
    arma::uvec sind = arma::stable_sort_index(simis);
    arma::uvec bstp = sind(bst);
    
    // most similar
    arma::mat  dbst = difng.rows(bstp);
    arma::vec  wgts = indst(bstp) / (sum(indst(bstp)));
    
    // prediction
    out.row(i) = f1.row(i);
    arma::rowvec delta = sum(dbst.each_col() % wgts, 0);
    out.row(i) += delta;
  }
  
  return out;
}