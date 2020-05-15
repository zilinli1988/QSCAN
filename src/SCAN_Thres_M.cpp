// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declare maxL1
arma::vec maxL1(int p, int Lmax, int Lmin, arma::mat x, arma::vec weights, arma::mat Cov, int times, int steplength);

// [[Rcpp::export]]
arma::vec SCAN_Thres_M(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, int steplength, arma::vec weights) {

	bool lower = false;
	bool logp = true;

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
	
	thr = maxL1(p, Lmax, Lmin, trans(x), weights, Cov, times, steplength);

	return thr;

}

arma::vec maxL1(int p, int Lmax, int Lmin, arma::mat x, arma::vec weights, arma::mat Cov, int times, int steplength)
{
	int ii, i, j, k, r;
	
	arma::vec sum0;
	sum0.zeros(times);
	
	arma::vec sumx;
	sumx.zeros(times);
	
	arma::vec threshold;
	threshold.zeros(times);
	
	bool lower = false;
	bool logp = true;
	
	double w = 0.0;

	int lengthnum = (Lmax-Lmin)/steplength + 1;
	

	for (ii = 0; ii < lengthnum; ii++)
	{
		i = Lmin + ii*steplength;
		
		for(r = 0; r < times; r++)
		{
			sum0[r] = 0;
		}

		w = 0;
		for(r = 0; r < times; r++)
		{
			for (k = 0; k < i; k++)
			{
				sum0[r] = sum0[r] + x(k,r)*weights(k);
			}
		}

		if(i>1)
		{
			w = arma::accu(Cov(arma::span(0, i - 1), arma::span(0, i - 1)));
		}
		if(i==1)
		{
			w = Cov(i - 1, i - 1);
		}
		
		for(r = 0; r < times; r++)
		{
			sumx[r] = pow(sum0[r],2)/w;
		}
		
		for(r = 0; r < times; r++)
		{
			if(sumx[r] > threshold[r])
			{
				threshold[r] = sumx[r];
			}
		}
		
		for (j = 1; j < (p - i + 1); j++)
		{
			for(r = 0; r < times; r++)
			{
				sum0[r] = sum0[r] - x(j - 1, r)*weights(j - 1) + x(j + i - 1, r)*weights(j + i - 1);
			}
			
			if (i > 1)
			{
				w = w - arma::accu(Cov(arma::span(j - 1, j - 1), arma::span(j, j + i - 2))) - arma::accu(Cov(arma::span(j, j + i - 2), arma::span(j - 1, j - 1))) - Cov(j - 1, j - 1);
				w = w + arma::accu(Cov(arma::span(j + i - 1, j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Cov(arma::span(j, j + i - 2), arma::span(j + i - 1, j + i - 1))) + Cov(j + i - 1, j + i - 1);
			}
			if (i == 1)
			{
				w = w - Cov(j - 1, j - 1) + Cov(j + i - 1, j + i - 1);
			}
			
			
			for(r = 0; r < times; r++)
			{
				sumx[r] = pow(sum0[r],2)/w;
			}
		
			for(r = 0; r < times; r++)
			{
				if(sumx[r]>threshold[r])
				{
					threshold[r] = sumx[r];
				}
			}
		}
	}

	return threshold;
}



