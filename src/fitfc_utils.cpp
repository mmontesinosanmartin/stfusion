// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat spatial_filtering(arma::mat& f1,
                            arma::mat& f2,
                            arma::ivec& fdims,
                            int& w,
                            arma::vec& wg,
                            int& nsim){
  
  int npxls = f2.n_rows;        // no. of pixels
  int clp = fdims[1] + (2 * w); // no. columns with the pad
  int cmx = fdims[1] + w;
  int rwp = fdims[0] + (2 * w); // no. rows with the pad
  int str = clp * w;
  int end = clp * (rwp - w);
  arma::uvec bsti = arma::linspace<arma::uvec>(0,nsim-1,nsim);
  arma::mat  f2ngb = arma::zeros(pow(2*w+1,2),f2.n_cols);
  arma::vec  dism = arma::zeros(pow(2*w+1,2));
  arma::mat  out(prod(fdims),f2.n_cols);
  
  // for each pixel
  for(int i = str; i < end; i++){
    // select the neighborhood
    int col = i  % clp;
    int row = floor(i /clp);
    
    if((col >= w) & (col < cmx)) {
      
      size_t z = 0;
      dism.fill(0);
      for(int j = -w; j <= w; j++) {   
        for(int k = -w; k <= w; k++) {
          int ind = (clp * j) + i + k;
          dism(z) = sum(pow(f1.row(ind) - f1.row(i), 2));
          f2ngb.row(z) = f2.row(ind);
          z++;
        }
      }
      
      // compute the weighted sum
      arma::uvec simpxl = arma::sort_index(dism);
      arma::uvec bstpxl = simpxl.elem(bsti);
      int outpi = (row-w) * fdims[1] + (col - w);
      out.row(outpi) = wg.elem(bstpxl).t()  / sum(wg.elem(bstpxl)) * f2ngb.rows(bstpxl);
    }
  }
  return(out);
}