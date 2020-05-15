// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declare maxL2
arma::vec maxL2(int p, int Lmax, int Lmin, arma::mat x, arma::vec weights, arma::mat Cov, int times);

// [[Rcpp::export]]
arma::vec Q_SCAN_Thres(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, arma::vec weights) {

	// variants number
	int p = G.n_cols;
	// sample size
	int n = G.n_rows;
	// covariates number
	int q = X.n_cols;
	// initial value of threshold
	arma::vec thr;
	thr.zeros(times);

	// pesudo-residual
	arma::mat y;
	y = arma::randn<arma::mat>(times,n);

	// pesudo-score
	arma::mat x;
	x.zeros(times,p);

	// t(X)*G
	arma::mat tX_G;
	tX_G.zeros(q,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);


	// Weights Matrix
	arma::mat W;
	W.zeros(p,p);


	if(fam == 0)
	{
		tX_G = trans(X)*G;
		Cov = trans(G)*G - trans(tX_G)*inv(trans(X)*X)*tX_G;
		x = (y*G - y*X*inv(trans(X)*X)*tX_G)*sigma;
	}else
	{
		tX_G = trans(X)*(arma::diagmat(working))*G;
		Cov = trans(G)*arma::diagmat(working)*G - trans(tX_G)*inv(trans(X)*arma::diagmat(working)*X)*tX_G;
		y = y*arma::diagmat(sqrt(working));
		x = (y*G - y*X*inv(trans(X)*arma::diagmat(working)*X)*tX_G)*sigma;
	}


	W.each_col() = weights;
	Cov = W%Cov;
	W.each_row() = trans(weights);
	Cov = Cov%W;
	Cov = Cov*pow(sigma,2);

	thr = maxL2(p, Lmax, Lmin, trans(x), weights, Cov, times);

	return thr;

}


