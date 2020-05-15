// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec maxL2(int p, int Lmax, int Lmin, arma::mat x, arma::vec weights, arma::mat Cov, int times)
{
	int ij, i, j, k, ii, jj, kk, r, s, ss;

	double w1 = 0.0;
	double w2 = 0.0;
	double c1 = 0.0;
	double c2 = 0.0;

	arma::vec sum0;
	sum0.zeros(times);

	arma::vec sumx;
	sumx.zeros(times);

	arma::vec threshold;
	threshold.zeros(times);

	for (ij = 0; ij < (Lmax - Lmin + 1); ij++)
	{
		i = Lmin + ij;

		for(r = 0; r < times; r++)
		{
			sum0[r] = 0;
		}

		w1 = 0.0;
		w2 = 0.0;
		c1 = 0.0;
		c2 = 0.0;

		for(r = 0; r < times; r++)
		{
			for (k = 0; k < i; k++)
			{
				sum0[r] = sum0[r] + pow(x(k,r),2)*pow(weights(k),2);
			}
		}

		for(k = 0; k < i; k++)
		{
			w1 = w1 + Cov(k, k);
		}

		for (ii = 0; ii < (i - 1); ii++)
		{
			for (jj = (ii + 1); jj < i; jj++)
			{
				w2 = w2 + Cov(ii, ii)*Cov(jj, jj) - Cov(ii, jj)*Cov(jj, ii);
			}
		}

		c1 = w1;
		c2 = pow(w1,2) - 2*w2;

		for(r = 0; r < times; r++)
		{
			sumx[r] = (sum0[r] - c1)/sqrt(2*c2);
			if(sumx[r] > threshold[r])
			{
				threshold[r] = sumx[r];
			}
		}

		for (j = 1; j < (p - i + 1); j++)
		{
			w1 = w1 - Cov(j - 1, j - 1) + Cov(j + i - 1, j + i - 1);
			for (kk = 1; kk < i; kk++)
			{
				w2 = w2 - (Cov(j - 1, j - 1)*Cov(j - 1 + kk, j - 1 + kk) - Cov(j - 1, j - 1 + kk)*Cov(j - 1 + kk, j - 1));
				w2 = w2 + (Cov(j + i - 1, j + i - 1)*Cov(j - 1 + kk, j - 1 + kk) - Cov(j + i - 1, j - 1 + kk)*Cov(j - 1 + kk, j + i - 1));
			}

			c1 = w1;
			c2 = pow(w1,2) - 2*w2;

			for(r = 0; r < times; r++)
			{
				sum0[r] = sum0[r] - pow(x(j - 1, r), 2)*pow(weights(j - 1), 2) + pow(x(j + i - 1, r), 2)*pow(weights(j + i - 1), 2);
				sumx[r] = (sum0[r]-c1)/sqrt(2*c2);

				if(sumx[r] > threshold[r])
				{
					threshold[r] = sumx[r];
				}

			}
		}
	}

	return threshold;
}



