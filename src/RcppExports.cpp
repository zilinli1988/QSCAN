// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Q_SCAN_Search
List Q_SCAN_Search(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec phenotype, arma::vec mu0, const double threshold, const int Lmax, const int Lmin, const int begid, const double f, arma::vec weights);
RcppExport SEXP _QSCAN_Q_SCAN_Search(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP phenotypeSEXP, SEXP mu0SEXP, SEXP thresholdSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP begidSEXP, SEXP fSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< const int >::type begid(begidSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_SCAN_Search(G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, begid, f, weights));
    return rcpp_result_gen;
END_RCPP
}
// Q_SCAN_Thres
arma::vec Q_SCAN_Thres(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, arma::vec weights);
RcppExport SEXP _QSCAN_Q_SCAN_Thres(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP timesSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_SCAN_Thres(G, X, working, sigma, fam, times, Lmax, Lmin, weights));
    return rcpp_result_gen;
END_RCPP
}
// SCAN_Search_M
List SCAN_Search_M(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec phenotype, arma::vec mu0, const double threshold, const int Lmax, const int Lmin, int steplength, arma::vec weights, const int begid, const double f);
RcppExport SEXP _QSCAN_SCAN_Search_M(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP phenotypeSEXP, SEXP mu0SEXP, SEXP thresholdSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weightsSEXP, SEXP begidSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type begid(begidSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(SCAN_Search_M(G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, steplength, weights, begid, f));
    return rcpp_result_gen;
END_RCPP
}
// SCAN_Thres_M
arma::vec SCAN_Thres_M(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, int steplength, arma::vec weights);
RcppExport SEXP _QSCAN_SCAN_Thres_M(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP timesSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(SCAN_Thres_M(G, X, working, sigma, fam, times, Lmax, Lmin, steplength, weights));
    return rcpp_result_gen;
END_RCPP
}
// maxL2
arma::vec maxL2(int p, int Lmax, int Lmin, arma::mat x, arma::vec weights, arma::mat Cov, int times);
RcppExport SEXP _QSCAN_maxL2(SEXP pSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP xSEXP, SEXP weightsSEXP, SEXP CovSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(maxL2(p, Lmax, Lmin, x, weights, Cov, times));
    return rcpp_result_gen;
END_RCPP
}
// regionfilter
arma::mat regionfilter(arma::mat candidate, const double f);
RcppExport SEXP _QSCAN_regionfilter(SEXP candidateSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type candidate(candidateSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(regionfilter(candidate, f));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_QSCAN_Q_SCAN_Search", (DL_FUNC) &_QSCAN_Q_SCAN_Search, 13},
    {"_QSCAN_Q_SCAN_Thres", (DL_FUNC) &_QSCAN_Q_SCAN_Thres, 9},
    {"_QSCAN_SCAN_Search_M", (DL_FUNC) &_QSCAN_SCAN_Search_M, 14},
    {"_QSCAN_SCAN_Thres_M", (DL_FUNC) &_QSCAN_SCAN_Thres_M, 10},
    {"_QSCAN_maxL2", (DL_FUNC) &_QSCAN_maxL2, 7},
    {"_QSCAN_regionfilter", (DL_FUNC) &_QSCAN_regionfilter, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_QSCAN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
