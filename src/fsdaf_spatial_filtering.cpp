// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "POOL.h"
using namespace Rcpp;

// [[Rcpp::export]]
double get_dist(int ind, int ref, int w,  int nrows, int ncols) {
  // reference position
  int row0 = floor(ref / ncols);
  int col0 = ref - (ncols * row0);
  // neighbor position
  int row1 = floor(ind / ncols);
  int col1 = ind - (ncols * row1);
  // distance
  double out = sqrt(pow(row0 - row1, 2) + pow(col0 - col1, 2));
  // return output
  return out;
}

// [[Rcpp::export]]
arma::mat fsdaf_spatial_filtering(arma::mat& c1,
                                  arma::mat& c2,
                                  arma::mat& f1,
                                  arma::mat& delta,
                                  arma::uvec dims,
                                  int w,
                                  int nsim,
                                  int nclass,
                                  double minr,
                                  double maxr){
  
  // initialize output
  int npx = c1.n_rows;
  int nbd = c1.n_cols;
  arma::mat out(npx, nbd);
  // thresholds
  arma::rowvec thr = arma::stddev(f1, 0, 0) * 2/nclass;
  // for each pixel
  for(int i = 0; i < npx; i++){
  // int i = 0;
    // stop if needed
    Rcpp::checkUserInterrupt();
    // get the neighbors
    arma::uvec ngbs = get_ngbs(i, w, dims[0], dims[1]);
    int nngbs = ngbs.n_elem;
    arma::mat c1ng = c1.rows(ngbs);
    arma::mat c2ng = c2.rows(ngbs);
    arma::mat f1ng = f1.rows(ngbs);
    // neighbor properties
    arma::vec simi(nngbs);
    arma::uvec inds(nngbs);
    arma::vec dist(nngbs);
    size_t z = 0;
    // for each neighbor
    for(int j  = 0; j < nngbs; j++){
      // compute similiarity (spectral distance)
      arma::rowvec sdis = abs(f1.row(ngbs[j]) - f1.row(i));
      arma::rowvec rdis = sdis/(f1.row(i) + 1e-4);
      // if below threshold
      if(all((sdis - thr.as_row()) < 0)){
        simi[z] = sum(rdis);
        inds[z] = ngbs[j];
        dist[z] = get_dist(i, ngbs[j], w, dims[0], dims[0]);
        z++;
      }
    }
    // leaving out unselected options
    simi = simi(arma::span(0, z - 1));
    inds = inds(arma::span(0, z - 1));
    dist = dist(arma::span(0, z - 1));
    // selection
    double spwgt = 10.;
    arma::uvec opts = {inds.n_elem, nsim};
    int ncand = opts.min();
    arma::uvec sind = arma::sort_index(simi % (dist/w + spwgt), "ascending");
    arma::uvec isel = inds(sind(arma::span(0,ncand - 1)));
    arma::vec dsel = dist(sind(arma::span(0,ncand - 1)));
    arma::vec ssel = simi(sind(arma::span(0,ncand - 1)));
    // // weights
    arma::vec iwg = (ssel + 1) % (1 + dsel/w);
    arma::vec wgt = (1/iwg) / sum(1/iwg);
    // prediction
    arma::mat deltam = delta.rows(isel);
    arma::mat deltaw = deltam.each_col() % wgt;
    arma::rowvec predic = f1.row(i) + sum(deltaw, 0);
    arma::rowvec naivep = f1.row(i) + (c2.row(i) - c1.row(i));
    out.row(i) = predic;
    for(int j = 0; j < nbd; j++) {
      if(out(i,j) < minr | out(i,j) > maxr) out(i,j) = naivep(j);
    }
  }
  // return result
  return out;
}
