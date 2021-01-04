#ifndef POOL_H
#define POOL_H

#include <RcppArmadillo.h>
arma::uvec get_ngbs(int i, int w, int nrow, int ncol);

#include <RcppArmadillo.h>
arma::uvec complete_obs(arma::mat m, int dimn);

#include <RcppArmadillo.h>
arma::uvec sequence(int n);

#include <RcppArmadillo.h>
int filter_sin(arma::mat x, arma::vec y);

#endif