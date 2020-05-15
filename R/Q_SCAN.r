#' Quadratic SCAN (Q-SCAN) statistic based procedure
#'
#' The \code{Q_SCAN} function takes in genotype, phenotype and covariates and detect the association between a
#' quantitative/dichotomous phenotype and a variant-set in a sequence by using quadratic scan statistic based procedure.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence, where n is the sample
#' size and p is the number of variants.
#' @param phenotype an n*1 phenotype vector, where n is the sample size.
#' @param X an n*q covariates matrix, where n is the sample size and q is the number of covariates.
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' Can be either "gaussian" for continuous phenotype or "binomial" for binary phenotype.
#' @param Lmax maximum number of variants in searching windows.
#' @param Lmin minimum number of variants in searching windows.
#' @param steplength difference of number of variants in searching windows, that is, the number of variants in
#' searching windows are Lmin, Lmin+steplength, Lmin+steplength,...,Lmax. (Default is 1).
#' @param times a number of pesudo-residuals (default = 2000).
#' @param alpha familty-wise/genome-wide significance level. (Default is 0.05).
#' @param f an overlap fraction, which controls for the overlapping proportion of of detected regions. For example,
#' when f=0, the detected regions are non-overlapped with each other,
#' and when f=1, we keep every susceptive region as detected regions. (Default is 0.)
#' @return The function returns a list with the following members:
#' @return \code{SCAN_res}:   A matrix that summarized the significant region detected by Q-SCAN.
#' The first column is the quadratic scan statistic of the detected region.
#' The next two columns are the location of the detected region (in sense of variants order).
#' The last column is the family-wise/genome-wide error rate of the detected region.
#' The result (0,0,0,1) means there is no significant region.
#' @return \code{SCAN_top1}:  A vector of length 4 which summarized the top 1 region detected by Q-SCAN.
#' The first element is the quadratic scan statistic of the detected region
#' The next two elements are the location of the detected region (in sense of variants order).
#' The last element is the family-wise/genome-wide p-value.
#' @return \code{SCAN_thres}: Empirical threshold of Q-SCAN for controlling the family-wise type I error at alpha level.
#' @return \code{SCAN_thres_boot}: A vector of Monte Carlo simulation sample for generating the empirical threshold. The 1-alpha quantile of this vector is
#' the empirical threshold.
#' @export


Q_SCAN <- function(genotype,phenotype,X,family,Lmax,Lmin,steplength=1,times=2000,alpha=0.05,f=0)
{
  maf <- colMeans(genotype)/2
  ## crop the sequence into sub sequence
  folds <- floor(dim(genotype)[2]/4000)

	samplesize <- dim(genotype)[1]
	# rank normal transformation
	if(family=="gaussian")
	{
		lmnull <- lm(phenotype~-1+X)
		sigma <- summary(lmnull)$sigma

		fam <- 0
		working <- rep(1,samplesize)

		mu0 <- lmnull$fitted
	}
	if(family!="gaussian")
	{
		# fit global null model
		glmnull <- glm(phenotype~-1+X,family=family)
		sigma <- sqrt(summary(glmnull)$dispersion)

		fam <- 1
		working <- glmnull$weights

		mu0 <- glmnull$fitted
	}


    ## SCANG-O
	L20 <- matrix(0,folds,times)
	res <- c()
	resmost <- c()

	weights <- dbeta(maf,1,1)
	subnum <- floor(dim(genotype)[2]/folds)

	for(i in 1:folds)
	{

		if(i<folds)
		{
			genotypesub <- genotype[,(subnum*(i-1)+1):(i*subnum+Lmax)]
			weightssub <- weights[(subnum*(i-1)+1):(i*subnum+Lmax)]
		}
		if(i==folds)
		{
			genotypesub <- genotype[,(subnum*(i-1)+1):dim(genotype)[2]]
			weightssub <- weights[(subnum*(i-1)+1):dim(genotype)[2]]
		}



		set.seed(19880615)
		threstemp <- Q_SCAN_Thres(genotypesub,X,working,sigma,fam,times,Lmax,Lmin,weightssub)

		begid <- subnum*(i-1)+1

		##### SCANG-O
		L20[i,] <- threstemp

		emL20 <- apply(L20,2,max)
		th0 <- quantile(emL20,1-alpha)

		restemp <- Q_SCAN_Search(genotypesub,X,working,sigma,fam,phenotype,mu0,th0,Lmax,Lmin,begid,f,weightssub)
		res <- rbind(res,restemp$res)
		resmost <- rbind(resmost,restemp$resmost)


	}

	rm(L20)
	gc()

	rm(mu0)
	gc()

	## SCANG-O
	res <- res[res[,1]>th0,]
	if(length(res)==0)
	{
		res <- c(0,0,0,1)
	}
	if(length(res)>4)
	{
		res <- regionfilter(res,f)
		if(length(res)==4)
		{
		  res[4] <- mean(emL20>res[1])
		}else
		{
		  res[,4] <- apply(res,1, function(z) mean(emL20>z[1]))
		}
	}

	mostnum <- which.max(resmost[,1])
	resmost <- resmost[mostnum,]
	resmost[4] <- mean(emL20>resmost[1])

	Lst <- list(SCAN_res=res,SCAN_top1=resmost,SCAN_thres=th0,SCAN_thres_boot=emL20)

	return(Lst)
}
