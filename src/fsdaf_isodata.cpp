// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

double modulus(int a, int n){
  return a - floor(a/n)*n;
}

arma::mat dist_to_class(arma::mat img, arma::mat clsc, int k){
  arma::mat out(img.n_rows, k);
  for(int i = 0; i < k; i++){
    out.col(i) = sqrt(sum(pow(img.each_row() - clsc.row(i),2), 1));
  }
  return out;
}

arma::mat class_stats(arma::mat img, arma::uvec m, int k){
  int nbd = img.n_cols;
  arma::mat out(k, (nbd + 1));
  for(int i = 0; i < k; i++){
    arma::uvec indx = find(m == i);
    out(i, arma::span(0,nbd - 1)) = sum(img.rows(indx), 0);
    out(i, nbd) = indx.n_elem;
  }
  return out;
}

// K = number of clusters
// I = number of iterations
// P = number of cluster pairs that can be merged
// minS = minimum number of points in a cluster
// maxStd = maximum value of standard deviation for splitting
// minDis = pairwise distance threshold to merge clusters
// M = threshold change in the cluster between each iteration.
// [[Rcpp::export]]
arma::uvec classify_isodata(arma::mat img,
                            int K = 4,
                            int I = 100,
                            int P = 2, 
                            double maxStdv = 500,
                            double minDis = 500,
                            double minS = 200,
                            double M = 0.05,
                            int k = 4) {
  
  // initialize output
  int npx = img.n_rows;
  int nbd = img.n_cols;
  arma::uvec m(npx);
  // class centers (initialize)
  arma::mat c = img.rows(arma::span(0,k - 1));
  
  // for each iteration
  for(int iter = 0; iter < I; iter++) {
    // first classification
    arma::mat distance(npx, k);
    arma::mat csum = arma::zeros(k, nbd);
    arma::vec cnum = arma::zeros(k);
    // distance to each class
    distance = dist_to_class(img, c, k);
    // minimum distance
    m = index_min(distance, 1);
    // class statistics
    arma::mat stats = class_stats(img, m, k);
    csum = stats.cols(arma::span(0, nbd - 1));
    cnum = stats.col(nbd);
    // handling infrequent class
    if(any(cnum < minS)){
      arma::vec delet = arma::zeros(k);
      int t = k;
      // check if few obs
      for(int i = 0; i < t; i++){
        // if so, delete
        if(cnum[i] < minS){
          k = k - 1;
          delet[i] = 1;
        }
      }
      arma::mat tdistance(npx, k);
      arma::mat tc(k, nbd);
      int ti = 0;
      // recalculate distance
      for(int i = 0; i < k; i++){
        if(delet[i] == 0){
          tdistance.col(ti) = distance.col(i);
          tc.row(ti) = c.row(i);
          ti += 1;
        }
      }
      distance = tdistance;
      c = tc;
      // recalculate clusters
      m = index_min(distance, 1);
      // recalculate statistics
      arma::mat stats = class_stats(img, m, k);
      csum = stats.cols(arma::span(0, nbd - 1));
      cnum = stats.col(nbd);
    }
    
    // every two, consider splitting classes
    if((modulus(iter, 2) == 0) || (k < K/2)){
      bool bsplit = false;
      double tmaxStd = -1;
      int nfeature = 0;
      int nclas = 0;
      // is std very high?
      for(int i = 0; i < k; i++){
        arma::uvec tindx = find(m == i);
        arma::mat  tvals = img.rows(tindx);
        arma::rowvec std = arma::stddev(tvals, 0);
        if(any(std > maxStdv)){
          int tnfeature = arma::index_max(std);
          if(std[tnfeature] > tmaxStd){
            tmaxStd = std[tnfeature];
            nfeature = tnfeature;
            nclas = i;
            bsplit = true;
          }
        }
      }
      // if so, split new cluster
      if(bsplit){
        // std of the class with highest variability
        arma::uvec splitindx = find(m == nclas);
        arma::mat  splitvals = img.rows(splitindx);
        arma::rowvec std = arma::stddev(splitvals, 0);
        // splitting into two
        arma::rowvec trow1 = c.row(nclas) - (M * std[nfeature]);
        arma::rowvec trow2 = c.row(nclas) + (M * std[nfeature]);
        k += 1;
        // replace the center
        c.row(nclas) = trow1;
        c = arma::join_cols(c, trow2);
        // update distances
        distance = dist_to_class(img, c, k);
        // update class
        m = index_min(distance, 1);
        // update stats
        arma::mat stats = class_stats(img, m, k);
        csum = stats.cols(arma::span(0, nbd - 1));
        cnum = stats.col(nbd);
      }
    }
    
    // every other , consider merging classes
    if((modulus(iter, 2) == 1) || (k > (K * 2))){
      bool bmerge = false;
      // similarity between clusters
      arma::mat eudist(k, k);
      eudist = dist_to_class(c, c, k);
      // a matrix defining which ones to merge
      int theight = round(k * (k - 1)/2);
      arma::mat distStruct(theight, 5);
      int cursor = 0;
      for(int classi = 0; classi < k; classi++){
        for(int classj = 0; classj < classi; classj++){
          distStruct.row(cursor) = {eudist(classi, classj), classi, classj, 0, 0};
          cursor += 1;
        }
      }
      // sort clusters based on similarity
      distStruct = distStruct.rows(arma::sort_index(distStruct.col(0), "ascend"));
      // if similarity between classes is below threshold, mark them
      for(int i = 0; i < theight; i++){
        if((distStruct(i, 4) == 0) && (distStruct(i,0) < minDis)){
          bmerge = true;
          distStruct(i,3) = 1;
          for(int j = 0; j < theight; j++){
            if((distStruct(j,1) == distStruct(i,1)) || (distStruct(j,2) == distStruct(i,1)) || (distStruct(j,1) == distStruct(i,2)) || (distStruct(j,2) == distStruct(i,2))){
              distStruct(j,4) = 1;
            }
          }
        }
      }
      arma::mat tc = c;
      // create a new merged class
      bool marker = false;
      for(int i = 0; i < theight; i++){
        if(distStruct(i, 3) == 1){
          int classa = (int) distStruct(i, 1);
          int classb = (int) distStruct(i, 2);
          tc.row(classa) = (tc.row(classa) + tc.row(classb))/2;
          tc.row(classb).fill(0);
          marker = true;
          k = k - 1;
        }
      }
      // if merged, recompute
      if(marker){
        arma::vec valid = arma::ones(tc.n_rows);
        for(int i = 0; i < tc.n_rows; i++){
          if(tc.row(i).is_zero()) valid[i] = 0;
        }
        c = tc.rows(find(valid == 1));
        // update distances, classes, and stats
        distance = dist_to_class(img, c, k);
        // minimum distance
        m = index_min(distance, 1);
        // class statistics
        arma::mat stats = class_stats(img, m, k);
        csum = stats.cols(arma::span(0, nbd - 1));
        cnum = stats.col(nbd);
      }
    }
    if(any(cnum == 0)){
      csum.replace(0,0.01);
      cnum.replace(0, 1);
    }
    arma::mat t = csum.each_col() /cnum;
    if(all(vectorise(t) == vectorise(c))){
      return m;
    }
    c = t;
  }
  // return classification
  return m;
}