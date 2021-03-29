// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// complete_obs
arma::uvec complete_obs(arma::mat m, int dimn);
RcppExport SEXP _stfusion_complete_obs(SEXP mSEXP, SEXP dimnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type dimn(dimnSEXP);
    rcpp_result_gen = Rcpp::wrap(complete_obs(m, dimn));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rmse
arma::rowvec cpp_rmse(arma::mat x, arma::mat y, bool byband);
RcppExport SEXP _stfusion_cpp_rmse(SEXP xSEXP, SEXP ySEXP, SEXP bybandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type byband(bybandSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rmse(x, y, byband));
    return rcpp_result_gen;
END_RCPP
}
// spatial_filtering
arma::mat spatial_filtering(arma::mat& f1, arma::mat& f2, arma::ivec& fdims, int& w, arma::vec& wg, int& nsim);
RcppExport SEXP _stfusion_spatial_filtering(SEXP f1SEXP, SEXP f2SEXP, SEXP fdimsSEXP, SEXP wSEXP, SEXP wgSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type fdims(fdimsSEXP);
    Rcpp::traits::input_parameter< int& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wg(wgSEXP);
    Rcpp::traits::input_parameter< int& >::type nsim(nsimSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_filtering(f1, f2, fdims, w, wg, nsim));
    return rcpp_result_gen;
END_RCPP
}
// get_ngbs
arma::uvec get_ngbs(int i, int w, int nrow, int ncol);
RcppExport SEXP _stfusion_get_ngbs(SEXP iSEXP, SEXP wSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ngbs(i, w, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// radio_par
arma::mat radio_par(arma::mat x, arma::mat y, arma::uvec dims, int w);
RcppExport SEXP _stfusion_radio_par(SEXP xSEXP, SEXP ySEXP, SEXP dimsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(radio_par(x, y, dims, w));
    return rcpp_result_gen;
END_RCPP
}
// starfm_spatial_filter
arma::mat starfm_spatial_filter(arma::mat& c1, arma::mat& c2, arma::mat& f1, arma::uvec dims, unsigned int redb, unsigned int nirb, int w, int ns);
RcppExport SEXP _stfusion_starfm_spatial_filter(SEXP c1SEXP, SEXP c2SEXP, SEXP f1SEXP, SEXP dimsSEXP, SEXP redbSEXP, SEXP nirbSEXP, SEXP wSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type redb(redbSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nirb(nirbSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(starfm_spatial_filter(c1, c2, f1, dims, redb, nirb, w, ns));
    return rcpp_result_gen;
END_RCPP
}
// apply_blur
arma::mat apply_blur(arma::mat& r, arma::uvec& rdims, double sigma);
RcppExport SEXP _stfusion_apply_blur(SEXP rSEXP, SEXP rdimsSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type rdims(rdimsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_blur(r, rdims, sigma));
    return rcpp_result_gen;
END_RCPP
}
// composite_genr
arma::vec composite_genr(arma::mat& x, arma::uvec& n);
RcppExport SEXP _stfusion_composite_genr(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(composite_genr(x, n));
    return rcpp_result_gen;
END_RCPP
}
// composite_lois
arma::uvec composite_lois(arma::mat& x, arma::vec& y, arma::ivec& cdims, int w);
RcppExport SEXP _stfusion_composite_lois(SEXP xSEXP, SEXP ySEXP, SEXP cdimsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type cdims(cdimsSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(composite_lois(x, y, cdims, w));
    return rcpp_result_gen;
END_RCPP
}
// filter_cor
int filter_cor(arma::mat x, arma::vec y);
RcppExport SEXP _stfusion_filter_cor(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(filter_cor(x, y));
    return rcpp_result_gen;
END_RCPP
}
// filter_sin
int filter_sin(arma::mat x, arma::vec y);
RcppExport SEXP _stfusion_filter_sin(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(filter_sin(x, y));
    return rcpp_result_gen;
END_RCPP
}
// local_lm
arma::mat local_lm(arma::mat& x, arma::mat& y, arma::uvec& n, arma::ivec& cdims, int w);
RcppExport SEXP _stfusion_local_lm(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP, SEXP cdimsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type cdims(cdimsSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(local_lm(x, y, n, cdims, w));
    return rcpp_result_gen;
END_RCPP
}
// spatial_filtering_par
arma::mat spatial_filtering_par(arma::mat ref, arma::mat img, arma::uvec dims, arma::uvec indx, int w, int nng, double spwgt);
RcppExport SEXP _stfusion_spatial_filtering_par(SEXP refSEXP, SEXP imgSEXP, SEXP dimsSEXP, SEXP indxSEXP, SEXP wSEXP, SEXP nngSEXP, SEXP spwgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ref(refSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type img(imgSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type nng(nngSEXP);
    Rcpp::traits::input_parameter< double >::type spwgt(spwgtSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_filtering_par(ref, img, dims, indx, w, nng, spwgt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stfusion_complete_obs", (DL_FUNC) &_stfusion_complete_obs, 2},
    {"_stfusion_cpp_rmse", (DL_FUNC) &_stfusion_cpp_rmse, 3},
    {"_stfusion_spatial_filtering", (DL_FUNC) &_stfusion_spatial_filtering, 6},
    {"_stfusion_get_ngbs", (DL_FUNC) &_stfusion_get_ngbs, 4},
    {"_stfusion_radio_par", (DL_FUNC) &_stfusion_radio_par, 4},
    {"_stfusion_starfm_spatial_filter", (DL_FUNC) &_stfusion_starfm_spatial_filter, 8},
    {"_stfusion_apply_blur", (DL_FUNC) &_stfusion_apply_blur, 3},
    {"_stfusion_composite_genr", (DL_FUNC) &_stfusion_composite_genr, 2},
    {"_stfusion_composite_lois", (DL_FUNC) &_stfusion_composite_lois, 4},
    {"_stfusion_filter_cor", (DL_FUNC) &_stfusion_filter_cor, 2},
    {"_stfusion_filter_sin", (DL_FUNC) &_stfusion_filter_sin, 2},
    {"_stfusion_local_lm", (DL_FUNC) &_stfusion_local_lm, 5},
    {"_stfusion_spatial_filtering_par", (DL_FUNC) &_stfusion_spatial_filtering_par, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_stfusion(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
