// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// GENERIC: functions useful for different purposes

// [[Rcpp::export]]
arma::uvec sequ(int n) {
  arma::uvec out = arma::linspace<arma::uvec>(0, n - 1, n);
  return out;
}

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
arma::vec ngb_dist(int i, arma::uvec inds, int nrow, int ncol){
  // reference row/col
  int row = floor(i / ncol);
  int col = i - (ncol * row);
  // pixel rows and columns
  arma::uvec urows = floor(inds / ncol);
  arma::uvec ucols = inds - (ncol * urows);
  arma::vec rows = arma::conv_to<arma::vec>::from(urows);
  arma::vec cols = arma::conv_to<arma::vec>::from(ucols);
  // distance
  arma::vec dist = sqrt(pow(rows - row,2) + pow(cols - col,2));
  return dist;
}


// CLOUD CLEANING: functions that help to remove cloudy pixels

// [[Rcpp::export]]
arma::vec distance_to_na(arma::vec& r, arma::ivec& dims, int maxdist){
  // initialize
  int npx = r.n_rows;
  arma::vec out(npx);
  out.fill(-1);
  // first fill
  arma::uvec nnas = arma::find_finite(r);
  out.rows(nnas).fill(maxdist);
  // parameters
  int nrows = dims(0);
  int ncols = dims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // if it is an edge
    if(r(i) == 1){
      // find raw neighbors
      arma::uvec ngbr = get_ngbs(i, maxdist, nrows, ncols);
      arma::vec  vngb = r.rows(ngbr);
      // find valid  values
      arma::uvec ngnna = arma::find_finite(vngb);
      arma::uvec ngbs = ngbr(ngnna);
      // compute distance
      arma::vec  ngds = ngb_dist(i, ngbs, nrows, ncols);
      arma::vec  outd = out(ngbs);
      arma::mat  dstc = join_rows(outd,ngds);
      // minimum distance
      out.rows(ngbs) = min(dstc, 1);
    }
  }
  return out; 
}



// UNSHARP MASKING: functions to de-blur coarse images

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
arma::mat blur(arma::mat& r, arma::uvec& rdims, double sigma) {
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

// [[Rcpp::export]]
arma::mat unsharp_masking(arma::mat& r, arma::uvec& rdims, double w) {
  // initialize output
  int npx = r.n_rows;
  int nbd = r.n_cols;
  arma::mat out(npx, nbd);
  // inner parameters
  int nrow = rdims(0);
  int ncol = rdims(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // int i = 0;
    // neighbors
    arma::uvec ngbs = get_ngbs(i, 1, nrow, ncol);
    arma::mat  rngb = r.rows(ngbs);
    // get the filter
    arma::uvec cnt = find(ngbs == i);
    rngb.shed_row(cnt[0]);
    // apply filter
    out.row(i) = r.row(i) + w * (r.row(i) - mean(rngb, 0));
  }
  return out;
}



// VIRTUAL IMAGE: generation through linear regression

// [[Rcpp::export]]
arma::uvec complete_obs(arma::mat m, int dimn, int hard){
  m.elem(find_finite(m)).fill(0);
  m.replace(arma::datum::nan, 1);
  arma::vec elems = sum(m, dimn).as_col();
  arma::uvec comp = find(elems <= hard);
  return comp;
}

// [[Rcpp::export]]
List filter_nna(arma::mat m){
  List out(2);
  // soft filter on rows
  // arma::uvec rfil = complete_obs(m, 1, m.n_cols - 2);
  // arma::mat  mrow = m.rows(rfil);
  // hard filter on columns
  // arma::mat  mxvl = mrow.cols(1, mrow.n_cols - 1);
  arma::mat  mxvl = m.cols(1, m.n_cols - 1);
  arma::uvec vfil = complete_obs(mxvl, 0, 0);
  // arma::mat  mfil = join_rows(mrow.col(0),mxvl.cols(vfil));
  arma::mat  mfil = join_rows(m.col(0),mxvl.cols(vfil));
  // return
  out[0] = vfil;
  out[1] = mfil;
  return out;
}

// [[Rcpp::export]]
int filter(arma::mat x, arma::vec y){
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
arma::mat local_regression(arma::mat& x, arma::mat& y, arma::ivec& cdims, int w) {
  // initialize output
  int npx = x.n_rows;
  int ntm = x.n_cols;
  int ncl = ntm + 1;
  arma::mat out(npx, ncl);
  // inner parameters
  int nrow = cdims(0);
  int ncol = cdims(1);
  arma::uvec kind(1);
  // for each pixel
  for(int i = 0; i < npx; i++) {
  // int i = 3834;
    // neighboring pixels
    arma::uvec inds = get_ngbs(i, w, nrow, ncol);
    // regression data
    arma::mat xng = x.rows(inds);
    arma::vec kns = arma::ones(xng.n_rows);
    arma::mat yng = y.rows(inds);
    // best time period
    int mxcor = filter(xng, yng);
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



// SPATIAL WEIGHTING: weighting predictions for similar pixels

// [[Rcpp::export]]
arma::mat sp_pred(arma::mat ref,
                   arma::mat img,
                   arma::uvec dims,
                   int w,
                   int nng) {
  
  // initialize
  int npx = ref.n_rows;
  int nbd = ref.n_cols;
  arma::mat out(npx, nbd);
  
  // parameters
  int nrows = dims(0);
  int ncols = dims(1);
  arma::uvec bst = sequ(nng);
  
  // for each pixel
  for(int i = 0; i < npx; i++) {
    // int i = 1699375;
    // stop if needed
    Rcpp::checkUserInterrupt();
    
    // current row-col
    int row = floor(i / ncols);
    int col = i - (ncols * row);
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
        int ind = (ncols * j) + i + k;
        simis(z) = sum(pow((ref.row(ind) - ref.row(i)), 2));
        indst(z) = 1. / (1 + sqrt(pow(j, 2) + pow(k, 2)) / w);
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
  // result
  return out;
}