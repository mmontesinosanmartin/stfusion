# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Blurs an image
#' 
#' @description blurs an image using a Gaussian kernel, which is a usual
#' way to describe the point spread function
#' 
#' @param r a matrix with pixel by bands dimensions
#' @param rdims a vector specifying the dimensions of the image
#' @param sigma standard deviation of the Gaussian kernel, usually equal to 1
#' 
#' @retruns a blurred image as a matrix
#' 
apply_blur <- function(r, rdims, sigma) {
    .Call(`_stfusion_apply_blur`, r, rdims, sigma)
}

#' @title Finds rows/columns without NAs
#' 
#' @description provides the index for the rows and columns without any
#' missing value
#' 
#' @param m a dataset, as a matrix
#' @param dimn can be 0 or 1 to find NAs across rows and columns respectively
#'  
complete_obs <- function(m, dimn) {
    .Call(`_stfusion_complete_obs`, m, dimn)
}

cpp_rmse <- function(x, y, byband) {
    .Call(`_stfusion_cpp_rmse`, x, y, byband)
}

spatial_filtering <- function(f1, f2, fdims, w, wg, nsim) {
    .Call(`_stfusion_spatial_filtering`, f1, f2, fdims, w, wg, nsim)
}

get_ngbs <- function(i, w, nrow, ncol) {
    .Call(`_stfusion_get_ngbs`, i, w, nrow, ncol)
}

filter_sin <- function(x, y) {
    .Call(`_stfusion_filter_sin`, x, y)
}

local_regression <- function(x, y, cdims, w) {
    .Call(`_stfusion_local_regression`, x, y, cdims, w)
}

corr_composite <- function(x, y, cdims, w) {
    .Call(`_stfusion_corr_composite`, x, y, cdims, w)
}

local_lm <- function(x, y, n, cdims, w) {
    .Call(`_stfusion_local_lm`, x, y, n, cdims, w)
}

composite <- function(x, n) {
    .Call(`_stfusion_composite`, x, n)
}

#' @title Calculate radiometric parameters
#' 
#' @description Calculates the slope and intercept to radiometrically correct
#' satellite images.
#' 
#' @details the function has to be applied for each band independently. 
#' Matrices x and y must have the pixels as rows and different dates as
#' columns. Both matrices must have the same length and NAs are not
#' allowed.
#' 
#' @param x a matrix with the original image (pixels x dates)
#' @param y a matrix with the reference image (pixels x dates)
#' @param dims a vector with the dimensions of the image (rows x columns)
#' @param w a number specifying the window radius (side 2w+1)
#' 
#' @returns a matrix with pixels x parameters (slope and intercept)
#' 
radio_par <- function(x, y, dims, w) {
    .Call(`_stfusion_radio_par`, x, y, dims, w)
}

sequence <- function(n) {
    .Call(`_stfusion_sequence`, n)
}

sp_pred <- function(ref, img, dims, w, nng, spwgt) {
    .Call(`_stfusion_sp_pred`, ref, img, dims, w, nng, spwgt)
}

vclamp <- function(x, minx, maxx) {
    .Call(`_stfusion_vclamp`, x, minx, maxx)
}

sp_pred_par <- function(ref, img, dims, indx, w, nng, spwgt) {
    .Call(`_stfusion_sp_pred_par`, ref, img, dims, indx, w, nng, spwgt)
}

