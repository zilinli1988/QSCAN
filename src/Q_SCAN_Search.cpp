// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Q_SCAN_Search(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec phenotype, arma::vec mu0, const double threshold, const int Lmax, const int Lmin, const int begid, const double f, arma::vec weights)
{
  int ij, i, j, k, ii, jj, kk, ss;
  int num = 0;
  double sum0 = 0.0;
  double w1 = 0.0;
  double w2 = 0.0;
  double c1 = 0.0;
  double c2 = 0.0;
  double sumx = 0.0;

  arma::rowvec x = trans(phenotype-mu0)*G;

  arma::mat candidate(1, 4);
  candidate.zeros();

  double summax = -100000.0;
  arma::mat candidatemax(1,4);
  candidatemax.zeros();

  int p = G.n_cols;
  int n = G.n_rows;
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


  for (ij = 0; ij < (Lmax-Lmin+1); ij++)
  {
    i = Lmin + ij;

    sum0 = 0;
    w1 = 0.0;
    w2 = 0.0;
    c1 = 0.0;
    c2 = 0.0;

    for (k = 0; k < i; k++)
    {
      sum0 = sum0 + pow(x(k), 2) *pow(weights(k), 2);
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
    sumx = (sum0-c1)/sqrt(2*c2);

    if (sumx>threshold)
    {
       num = num + 1;
       candidate.resize(num, 4);
       candidate(num - 1, 0) = sumx;
       candidate(num - 1, 1) = 1;
       candidate(num - 1, 2) = i;
    }

    if(sumx>summax)
    {
       summax=sumx;
       candidatemax(0,0)=sumx;
       candidatemax(0,1)=1+begid-1;
       candidatemax(0,2)=i+begid-1;
    }


    for (j = 1; j < (p - i + 1); j++)
    {
      sum0 = sum0 - pow(x(j - 1), 2)*pow(weights(j - 1), 2) + pow(x(j + i - 1), 2)*pow(weights(j + i - 1), 2);
      w1 = w1 - Cov(j - 1, j - 1) + Cov(j + i - 1, j + i - 1);
      for (kk = 1; kk < i; kk++)
      {
        w2 = w2 - (Cov(j - 1, j - 1)*Cov(j - 1 + kk, j - 1 + kk) - Cov(j - 1, j - 1 + kk)*Cov(j - 1 + kk, j - 1));
        w2 = w2 + (Cov(j + i - 1, j + i - 1)*Cov(j - 1 + kk, j - 1 + kk) - Cov(j + i - 1, j - 1 + kk)*Cov(j - 1 + kk, j + i - 1));
      }
	  
	  c1 = w1;
	  c2 = pow(w1,2) - 2*w2;
	  sumx = (sum0-c1)/sqrt(2*c2);


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
        candidatemax(0,1)=j + 1 + begid - 1;
        candidatemax(0,2)=j + i + begid - 1;
      }
    }
  }

  arma::uvec indices = sort_index(-candidate.col(0));
  candidate = candidate.rows(indices);

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
