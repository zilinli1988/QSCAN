# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Q_SCAN_Search <- function(G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, begid, f, weights) {
    .Call(`_QSCAN_Q_SCAN_Search`, G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, begid, f, weights)
}

Q_SCAN_Thres <- function(G, X, working, sigma, fam, times, Lmax, Lmin, weights) {
    .Call(`_QSCAN_Q_SCAN_Thres`, G, X, working, sigma, fam, times, Lmax, Lmin, weights)
}

SCAN_Search_M <- function(G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, steplength, weights, begid, f) {
    .Call(`_QSCAN_SCAN_Search_M`, G, X, working, sigma, fam, phenotype, mu0, threshold, Lmax, Lmin, steplength, weights, begid, f)
}

SCAN_Thres_M <- function(G, X, working, sigma, fam, times, Lmax, Lmin, steplength, weights) {
    .Call(`_QSCAN_SCAN_Thres_M`, G, X, working, sigma, fam, times, Lmax, Lmin, steplength, weights)
}

maxL2 <- function(p, Lmax, Lmin, x, weights, Cov, times) {
    .Call(`_QSCAN_maxL2`, p, Lmax, Lmin, x, weights, Cov, times)
}

regionfilter <- function(candidate, f) {
    .Call(`_QSCAN_regionfilter`, candidate, f)
}

