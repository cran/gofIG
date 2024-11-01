#' The first Allison-Betsch-Ebner-Visagie goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family due to Allison et al. (2019). Two different estimation procedures are implemented, namely the method of moment and the maximum likelihood method.
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param meth  method of estimation used. Possible values are \code{'MME'} for moment estimation and \code{'MLE'} for maximum likelihood estimation.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$parameter}}{the value of the tuning parameter.}
#'         \item{\code{$est.method}}{the estimation method used.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the inverse Gaussian distribution. The p value is obtained by a parametric bootstrap procedure.
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2019) "New weighted \eqn{L^2}-type tests for the inverse Gaussian distribution", arXiv:1910.14119. \href{https://arxiv.org/abs/1910.14119}{LINK}
#'
#' @examples
#' test.ABEV1(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.ABEV1 <- function(data,a=10,meth='MME',B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
	n   = length(data)
	muH = mean(data)
	if (meth=='MLE'){
	  lambdaH = 1/mean(1/data-1/muH)
	} else if (meth=='MME') {
	  lambdaH = muH^3/(mean(data^2)-muH^2)
	}
	phiH = lambdaH/muH
	TnH0 = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = ABEV1(Xstar,a,meth)
	}
	TnH0 = sort(TnH0)
	TnX  = ABEV1(data,a,meth)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "ABEV1","parameter" = a, "est.method"=meth, "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)
}



#' The second Allison-Betsch-Ebner-Visagie goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family due to Allison et al. (2019). Two different estimation procedures are implemented, namely the method of moment and the maximum likelihood method.
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param meth  method of estimation used. Possible values are \code{'MME'} for moment estimation and \code{'MLE'} for maximum likelihood estimation.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$parameter}}{the value of the tuning parameter.}
#'         \item{\code{$est.method}}{the estimation method used.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the inverse Gaussian distribution. The p value is obtained by a parametric bootstrap procedure.
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2019) "New weighted \eqn{L^2}-type tests for the inverse Gaussian distribution", arXiv:1910.14119. \href{https://arxiv.org/abs/1910.14119}{LINK}
#'
#' @examples
#' test.ABEV2(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.ABEV2 <- function(data,a=10,meth='MME',B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
  n   = length(data)
  muH = mean(data)
  if (meth=='MLE'){
    lambdaH = 1/mean(1/data-1/muH)
  } else if (meth=='MME') {
    lambdaH = muH^3/(mean(data^2)-muH^2)
  }
  phiH = lambdaH/muH
  TnH0 = rep(0,B)
  for (j in 1:B){
    Xstar   = rmutil::rinvgauss(n,1,1/phiH)
    TnH0[j] = ABEV2(Xstar,a,meth)
  }
  TnH0 = sort(TnH0)
  TnX  = ABEV2(data,a,meth)
  pval = mean(TnH0>TnX)
  result <-list("Test" = "ABEV2","parameter" = a, "est.method"=meth, "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
  attr(result, "class") <- "gofIG"
  return(result)
}



#' The first Henze-Klar goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$parameter}}{the value of the tuning parameter.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test statistics is a weighted integral over the squared modulus of some measure of deviation of the empirical distribution of given data from the family of inverse Gaussian laws, expressed by means of the empirical Laplace transform.
#'
#' @references
#' Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
#'
#' @examples
#' test.HK1(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.HK1 <- function(data,a=0,B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
	n       = length(data)
	muH     = mean(data)
	lambdaH = 1/mean(1/data-1/muH)
	phiH    = lambdaH/muH
	TnH0    = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = HK1(Xstar,a)
	}
	TnH0 = sort(TnH0)
	TnX  = HK1(data,a)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "HK1","parameter" = a, "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)
}



#' The second Henze-Klar goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
#'
#' @param data  a vector of positive numbers.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test statistic is a weighted integral over the squared modulus of some measure of deviation of the empirical distribution of given data from the family of inverse Gaussian laws, expressed by means of the empirical Laplace transform.
#'
#' @references
#' Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
#'
#' @examples
#' test.HK2(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.HK2 <- function(data,B){
  if (!is.vector(data)){warning("data must be of type vector")}
	n       = length(data)
	muH     = mean(data)
	lambdaH = 1/mean(1/data-1/muH)
	phiH    = lambdaH/muH
	TnH0    = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = HK2(Xstar)
	}
	TnH0 = sort(TnH0)
	TnX  = HK2(data)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "HK2", "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)
}



#' The Kolmogorov-Smirnov goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family in the spirit of Kolmogorov and Smirnov. Note that this tests the composite hypothesis of fit to the family of inverse Gaussian distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#' @details
#' The Kolmogorov Smirnov test is computed as described in Allison et. al. (2019). The p value is obtained by a parametric bootstrap procedure.
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2019) "New weighted \eqn{L^2}-type tests for the inverse Gaussian distribution", arXiv:1910.14119. \href{https://arxiv.org/abs/1910.14119}{LINK}
#'
#' @examples
#' test.KS(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.KS <- function(data,B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
	n       = length(data)
	muH     = mean(data)
	lambdaH = 1/mean(1/data-1/muH)
	phiH    = lambdaH/muH
	TnH0    = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = KS(Xstar)
	}
	TnH0 = sort(TnH0)
	TnX  = KS(data)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "KS", "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)
}



#' The Cramer-von Mises goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family in the spirit of Cramer and von Mises. Note that this tests the composite hypothesis of fit to the family of inverse Gaussian distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#' @details
#' The Cramer-von Mises test is computed as described in Allison et. al. (2019). The p value is obtained by a parametric bootstrap procedure.
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2019) "New weighted \eqn{L^2}-type tests for the inverse Gaussian distribution", arXiv:1910.14119. \href{https://arxiv.org/abs/1910.14119}{LINK}
#'
#' @examples
#' test.CM(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.CM <- function(data,B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
	n       = length(data)
	muH     = mean(data)
	lambdaH = 1/mean(1/data-1/muH)
	phiH    = lambdaH/muH
	TnH0    = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = CM(Xstar)
	}
	TnH0 = sort(TnH0)
	TnX  = CM(data)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "CM", "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)

}



#' The Anderson-Darling goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family in the spirit of Anderson and Darling. Note that this tests the composite hypothesis of fit to the family of inverse Gaussian distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#' @details
#' The Anderson-Darling test is computed as described in Allison et. al. (2019). The p value is obtained by a parametric bootstrap procedure.
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2019) "New weighted \eqn{L^2}-type tests for the inverse Gaussian distribution", arXiv:1910.14119. \href{https://arxiv.org/abs/1910.14119}{LINK}
#'
#' @examples
#' test.AD(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.AD <- function(data,B=500){
  if (!is.vector(data)){warning("data must be of type vector")}
	n       = length(data)
	muH     = mean(data)
	lambdaH = 1/mean(1/data-1/muH)
	phiH    = lambdaH/muH
	TnH0    = rep(0,B)
	for (j in 1:B){
		Xstar   = rmutil::rinvgauss(n,1,1/phiH)
		TnH0[j] = AD(Xstar)
	}
	TnH0 = sort(TnH0)
	TnX  = AD(data)
	pval = mean(TnH0>TnX)
	result <-list("Test" = "AD", "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
	attr(result, "class") <- "gofIG"
	return(result)
}



#' The Baringhaus-Gaigall goodness-of-fit test for the inverse Gaussian family
#'
#' @description
#' This function computes the goodness-of-fit test for the inverse Gaussian family due to Baringhaus and Gaigall (2015).
#'
#' @param data  a vector of positive numbers.
#' @param B     number of bootstrap iterations used to obtain p value.
#'
#' @return a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
#' \describe{
#'         \item{\code{$Test}}{the name of the used test statistic.}
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$p.value}}{the approximated p value.}
#'         \item{\code{$par.est}}{the estimated parameters.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @references
#' Baringhaus, L.  Gaigall, D. (2015). "On an independence test approach to the goodness-of-fit problem", Journal of Multivariate Analysis, 140, 193-208. \doi{https://doi.org/10.1016/j.jmva.2015.05.013}
#'
#' @examples
#' test.BG(rmutil::rinvgauss(20,2,1),B=100)
#'
#' @export
test.BG <- function(data,B){
  if (!is.vector(data)){warning("data must be of type vector")}
  n       = length(data)
  muH     = mean(data)
  lambdaH = 1/mean(1/data-1/muH)
  phiH    = lambdaH/muH
  TnH0    = rep(0,B)
  for (j in 1:B){
    Xstar   = rmutil::rinvgauss(n,1,1/phiH)
    TnH0[j] = BG(Xstar)
  }
  TnH0 = sort(TnH0)
  TnX  = BG(data)
  pval = mean(TnH0>TnX)
  result <-list("Test" = "BG", "T.value"=TnX,"p.value"=pval,"par.est"=c(muH,lambdaH),"boot.run"=B)
  attr(result, "class") <- "gofIG"
  return(result)
}
