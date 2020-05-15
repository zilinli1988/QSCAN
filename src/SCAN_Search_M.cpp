
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SCAN_Search_M(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec phenotype, arma::vec mu0, const double threshold, const int Lmax, const int Lmin, int steplength, arma::vec weights, const int begid,  const double f)
{

	int ij, i, j, k;
	int num = 0;

	double sum0 = 0.0;
	double sumw = 0.0;
	double sumx = 0.0;
	arma::mat candidate(1, 4);
	candidate.zeros();

	double summax = -100000.0;
	arma::mat candidatemax(1,4);
	candidatemax.zeros();

	arma::rowvec x = trans(phenotype-mu0)*G;

	int p = G.n_cols;
	// int n = G.n_rows;
	int q = X.n_cols;

	// t(X)*G
	arma::mat tX_G;
	tX_G.zeros(q,p);

	// Weights Matrix
	arma::mat W;
	W.zeros(p,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	if(fam == 0)
	{
		tX_G = trans(X)*G;
		Cov = trans(G)*G - trans(tX_G)*inv(trans(X)*X)*tX_G;
	}else
	{
		tX_G = trans(X)*(arma::diagmat(working))*G;
		Cov = trans(G)*arma::diagmat(working)*G - trans(tX_G)*inv(trans(X)*arma::diagmat(working)*X)*tX_G;
	}

	W.each_col() = weights;
	Cov = W%Cov;
	W.each_row() = trans(weights);
	Cov = Cov%W;
	Cov = Cov*pow(sigma,2);

	int lengthnum = (Lmax-Lmin)/steplength + 1;

	for (ij = 0; ij < lengthnum; ij++)
	{
		i = Lmin + ij*steplength;

		sum0 = 0;
		sumw = 0;
		for (k = 0; k < i; k++)
		{
			sum0 = sum0 + x(k)*weights(k);
		}
		if (i>1)
		{
			sumw = arma::accu(Cov(arma::span(0, i - 1), arma::span(0, i - 1)));
		}
		if (i == 1)
		{
			sumw = Cov(i - 1, i - 1);
		}

		sumx = pow(sum0, 2) / sumw;

		if (sumx > threshold)
		{
			num = num + 1;
			candidate.resize(num, 4);
			candidate(num - 1, 0) = sumx;
			candidate(num - 1, 1) = 1;
			candidate(num - 1, 2) = i;
		}
		if(sumx > summax)
		{
			summax = sumx;
			candidatemax(0,0) = sumx;
			candidatemax(0,1) = 1 + begid - 1;
			candidatemax(0,2) = i + begid - 1;
		}

		for (j = 1; j < (p - i + 1); j++)
		{
			sum0 = sum0 - x(j - 1)*weights(j - 1) + x(j + i - 1)*weights(j + i - 1);
			if (i > 1)
			{
				sumw = sumw - arma::accu(Cov(arma::span(j - 1, j - 1), arma::span(j, j + i - 2))) - arma::accu(Cov(arma::span(j, j + i - 2), arma::span(j - 1, j - 1))) - Cov(j - 1, j - 1);
				sumw = sumw + arma::accu(Cov(arma::span(j + i - 1, j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Cov(arma::span(j, j + i - 2), arma::span(j + i - 1, j + i - 1))) + Cov(j + i - 1, j + i - 1);
			}
			if (i == 1)
			{
				sumw = sumw - Cov(j - 1, j - 1) + Cov(j + i - 1, j + i - 1);
			}

			sumx = pow(sum0, 2) / sumw;

			if (sumx>threshold)
			{
				num = num + 1;
				candidate.resize(num, 4);
				candidate(num - 1, 0) = sumx;
				candidate(num - 1, 1) = j + 1;
				candidate(num - 1, 2) = j + i;
			}
			if(sumx>summax)
			{
				summax=sumx;
				candidatemax(0,0)=sumx;
				candidatemax(0,1)=j+1+begid-1;
				candidatemax(0,2)=j+i+begid-1;
			}

		}
	}

	arma::uvec indices = sort_index(-candidate.col(0));
	candidate = candidate.rows(indices);

	int ii,jj,kk;

	double loc_left = 0;
	double loc_right = 0;

	for (ii = 0; ii < (num-1); ii++)
	{
		if (candidate(ii,3) < 1)
		{
			for (jj = ii + 1; jj < num; jj++)
			{
				if(candidate(ii, 1) < candidate(jj, 1))
				{
					loc_left = candidate(jj, 1);
				}else
				{
					loc_left = candidate(ii, 1);
				}
				if(candidate(ii, 2) < candidate(jj, 2))
				{
					loc_right = candidate(ii, 2);
				}else
				{
					loc_right = candidate(jj, 2);
				}

				if (loc_right > loc_left - 1)
				{
					if((loc_right-loc_left + 1)/(candidate(jj,2) - candidate(jj,1) + 1) > f)
					{
						candidate(jj, 3) = 1;
					}


				}
			}
		}
	}

	int num1 = 0;
	arma::mat res(1, 4);
	res.zeros();
	for (kk = 0; kk<num; kk++)
	{
		if (candidate(kk, 3)<1)
		{
			num1 = num1 + 1;
			res.resize(num1, 4);
			res(num1 - 1, 0) = candidate(kk, 0);
			res(num1 - 1, 1) = candidate(kk, 1) + begid - 1;
			res(num1 - 1, 2) = candidate(kk, 2) + begid - 1;
		}
	}
	return List::create(Named("res") = res, Named("resmost") = candidatemax);

}

