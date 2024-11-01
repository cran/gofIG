#' The first Allison-Betsch-Ebner-Visagie test statistic
#'
#' @description
#' This function computes the first test statistic of the goodness-of-fit tests for the inverse Gaussian family due to Allison et al. (2022). Two different estimation procedures are implemented, namely the method of moment and the maximum likelihood method.
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param meth  method of estimation used. Possible values are \code{'MME'} for moment estimation and \code{'MLE'} for maximum likelihood estimation.
#'
#' @return value of the test statistic.
#' 
#'
#' @details 
#' The numerically stable test statistic for the first Allison-Betsch-Ebner-Visagie test is defined as:
#' \deqn{ABEV1_{n,a} = \frac{1}{4n} \sum_{j,k=1}^{n} \left( \hat{\varphi}_n + \frac{3}{Y_{n,j}} - \frac{\hat{\varphi}_n}{Y_{n,j}^2} \right) \left( \hat{\varphi}_n + \frac{3}{Y_{n,k}} - \frac{\hat{\varphi}_n}{Y_{n,k}^2} \right) h_{1,a}(Y_{n,j}, Y_{n,k})}
#' \deqn{- 2 \left( \hat{\varphi}_n + \frac{3}{Y_{n,j}} - \frac{\hat{\varphi}_n}{Y_{n,j}^2} \right) h_{2,a}(Y_{n,j}, Y_{n,k})}
#' \deqn{- 2 \left( \hat{\varphi}_n + \frac{3}{Y_{n,k}} - \frac{\hat{\varphi}_n}{Y_{n,k}^2} \right) h_{2,a}(Y_{n,k}, Y_{n,j})}
#' \deqn{+ \frac{4}{a} e^{-a \max(Y_{n,j}, Y_{n,k})},}
#' with \eqn{\hat{\varphi}_n = \frac{\hat{\lambda}_n}{\hat{\mu}_n}}, where \eqn{\hat{\mu}_n,\hat{\lambda}_n} are consistent estimators of \eqn{\mu, \lambda}, respectively, the parameters of the inverse Gaussian distribution. Furthermore \eqn{Y_{n,j} = \frac{X_j}{\hat{\mu}_n}}, \eqn{j = 1,...,n}, for \eqn{(X_j)_{j = 1,...,n}}, a sequence of  independent observations of a positive random variable \eqn{X}. The functions \eqn{h_{i,a}(s,t)}, \eqn{i = 1,2}, are defined in Allison et al. (2022), section 5.1.  
#' The null hypothesis is rejected for large values of the test statistic \eqn{ABEV1_{n,a}}.
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2022) "On Testing the Adequacy of the Inverse Gaussian Distribution". \href{https://www.mdpi.com/2227-7390/10/3/350}{LINK}
#'
#' @examples
#' ABEV1(rmutil::rinvgauss(20,2,1),a=10,meth='MLE')
#' 
#' @export
ABEV1 <- function(data,a=10,meth='MME'){
  n   = length(data)
  muH = mean(data)
  if (meth=='MLE'){
    lambdaH = 1/mean(1/data-1/muH)
  } else if (meth=='MME') {
    lambdaH = muH^3/(mean(data^2)-muH^2)
  }
  phiH = lambdaH/muH
  Y    = data/muH
  sum1 = 0
  sum2 = 0
  sum3 = 0
  sum4 = 0
  for (j in 1:n){
    for (k in 1:n){
      sum1 = sum1+(phiH+3/Y[j]-phiH/(Y[j])^2)*(phiH+3/Y[k]-phiH/(Y[k])^2)*h1(Y[j],Y[k],a)
      sum2 = sum2+(phiH+3/Y[j]-phiH/(Y[j])^2)*h2(Y[j],Y[k],a)
      sum3 = sum3+(phiH+3/Y[k]-phiH/(Y[k])^2)*h2(Y[k],Y[j],a)
      sum4 = sum4+exp(-a*max(Y[j],Y[k]))/a
    }
  }
  TestStat = (1/(4*n))*(sum1-2*(sum2+sum3)+4*sum4)
  return(TestStat)
}



#' The second Allison-Betsch-Ebner-Visagie test statistic
#'
#' @description
#' This function computes the second test statistic of the goodness-of-fit tests for the inverse Gaussian family due to Allison et al. (2022). Two different estimation procedures are implemented, namely the method of moment and the maximum likelihood method.
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param meth  method of estimation used. Possible values are \code{'MME'} for moment estimation and \code{'MLE'} for maximum likelihood estimation.
#'
#' @return value of the test statistic.
#' 
#' @details 
#' The numerically stable test statistic for the second Allison-Betsch-Ebner-Visagie test is defined as:
#' \deqn{ABEV2_{n,a} = \frac{1}{4n} \sum_{j,k=1}^{n} \left( \hat{\varphi}_n + \frac{3}{Y_{n,j}} - \frac{\hat{\varphi}_n}{Y_{n,j}^2} \right) \left( \hat{\varphi}_n + \frac{3}{Y_{n,k}} - \frac{\hat{\varphi}_n}{Y_{n,k}^2} \right) \tilde{h}_{1,a}(Y_{n,j}, Y_{n,k})}
#' \deqn{- 2 \left( \hat{\varphi}_n + \frac{3}{Y_{n,j}} - \frac{\hat{\varphi}_n}{Y_{n,j}^2} \right) \tilde{h}_{2,a}(Y_{n,j}, Y_{n,k})}
#' \deqn{- 2 \left( \hat{\varphi}_n + \frac{3}{Y_{n,k}} - \frac{\hat{\varphi}_n}{Y_{n,k}^2} \right) \tilde{h}_{2,a}(Y_{n,k}, Y_{n,j})}
#' \deqn{+ 4 \frac{\sqrt{\pi}}{a} \Phi \left( - \sqrt{2a} \max(Y_{n,j}, Y_{n,k}) \right),}
#' with \eqn{\hat{\varphi}_n = \frac{\hat{\lambda}_n}{\hat{\mu}_n}}, where \eqn{\hat{\mu}_n,\hat{\lambda}_n} are consistent estimators of \eqn{\mu, \lambda}, respectively, the parameters of the inverse Gaussian distribution. Furthermore \eqn{Y_{n,j} = \frac{X_j}{\hat{\mu}_n}}, \eqn{j = 1,...,n}, for \eqn{(X_j)_{j = 1,...,n}}, a sequence of  independent observations of a positive random variable \eqn{X}. 
#' The functions \eqn{\tilde{h}_{i,a}(s,t)}, \eqn{i = 1,2}, are defined in Allison et al. (2022), section 5.1, and \eqn{\Phi} denotes the distribution function of the standard normal distribution.
#' The null hypothesis is rejected for large values of the test statistic \eqn{ABEV2_{n,a}}.
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2022) "On Testing the Adequacy of the Inverse Gaussian Distribution". \href{https://www.mdpi.com/2227-7390/10/3/350}{LINK}
#' 
#' @examples
#' ABEV2(rmutil::rinvgauss(20,2,1),a=10,meth='MLE')
#' 
#' @export
ABEV2 <- function(data,a=10,meth='MME'){
  n   = length(data)
  muH = mean(data)
  if (meth=='MLE'){
    lambdaH = 1/mean(1/data-1/muH)
  } else if (meth=='MME') {
    lambdaH = muH^3/(mean(data^2)-muH^2)
  }
  phiH = lambdaH/muH
  Y    = data/muH
  sum1 = 0
  sum2 = 0
  sum3 = 0
  sum4 = 0
  for (j in 1:n){
    for (k in 1:n){
      sum1 = sum1+(phiH+3/Y[j]-phiH/(Y[j])^2)*(phiH+3/Y[k]-phiH/(Y[k])^2)*h1t(Y[j],Y[k],a)
      sum2 = sum2+(phiH+3/Y[j]-phiH/(Y[j])^2)*h2t(Y[j],Y[k],a)
      sum3 = sum3+(phiH+3/Y[k]-phiH/(Y[k])^2)*h2t(Y[k],Y[j],a)
      sum4 = sum4+sqrt(pi/a)*stats::pnorm(-sqrt(2*a)*max(Y[j],Y[k]))
    }
  }
  TestStat = (1/(4*n))*(sum1-2*(sum2+sum3)+4*sum4)
  return(TestStat)
}



#' The first Henze-Klar test statistic
#'
#' @description
#' This function computes the first test statistic of the goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#'
#' @return value of the test statistic
#' 
#' @details 
#' The representation of the first Henze-Klar test statistic used for computation is given by: 
#' \deqn{HK_{n,a}^{(1)}= \frac{\hat{\varphi}_n}{n} \sum_{j,k=1}^{n} \hat{Z}_{jk}^{-1} \left\{ 1 - (Y_j + Y_k) \left( 1 + \sqrt{\frac{\pi}{2\hat{Z}_{jk}}} \text{erfce}\left( \sqrt{\frac{\hat{Z}_{jk}}{2}} \right) \right) + \left( 1 + \frac{2}{\hat{Z}_{jk}} \right) Y_j Y_k \right\},}
#' with \eqn{\hat{\varphi}_n = \frac{\hat{\lambda}_n}{\hat{\mu}_n}}, where \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution. 
#' Furthermore \eqn{\hat{Z}_{jk} = \hat{\varphi}_n(Y_j + Y_k +a)}, where \eqn{Y_i = \frac{X_i}{\hat{\mu}_n}} for \eqn{(X_i)_{i = 1,...,n}}, a sequence of  independent observations of a nonnegative random variable \eqn{X}.
#' To ensure numerical stability of the implementation the exponentially scaled complementary error function \eqn{\text{erfce}(x)} is used: \eqn{\text{erfce}(x) = \exp{(x^2)}\text{erfc}(x)}, with \eqn{\text{erfc}(x) = 2\int_x^\infty \exp{(-t^2)}dt/\pi}.
#' The null hypothesis is rejected for large values of the test statistic \eqn{HK_{n,a}^{(1)}}.
#' @references
#' Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
#'
#' @examples
#' HK1(rmutil::rinvgauss(20,2,1))
#' 
#' @export
HK1 <- function(data,a=0){
  n       = length(data)
  muH     = mean(data)
  lambdaH = 1/mean(1/data-1/muH)
  phiH    = lambdaH/muH
  Y       = data/muH
  Tmat    = matrix(0,n,n)
  for (j in 1:n){
    for (k in 1:n){
      ZjkH      = phiH*(Y[j]+Y[k]+a)
      T1a       = Y[j]+Y[k]
      T1b       = 1+sqrt(pi/(ZjkH*2))*erfce(sqrt(ZjkH/2))
      T2        = (1+2/ZjkH)*Y[j]*Y[k]
      Tmat[j,k] = (1 - T1a*T1b + T2)/ZjkH
    }
  }
  return(phiH/n*sum(Tmat))
}



#' The second Henze-Klar test statistic
#'
#' @description
#' This function computes the test statistic of  the second goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
#'
#' @param data  a vector of positive numbers.
#'
#' @return value of the test statistic.
#'
#' @details 
#' The representation of the second Henze-Klar test statistic used for computation \eqn{(a = 0)} is given by:
#' \deqn{HK_{n,0}^{(2)} = \frac{1}{n} \sum_{j,k=1}^{n} Z_{jk}^{-1} - 2 \sum_{j=1}^{n} Z_j^{-1} \left\{ 1 - \sqrt{\frac{\pi \hat{\varphi}_n}{2 Z_j}} \, \mathrm{erfce} \left( \frac{\hat{\varphi}_n^{1/2} (Z_j + 1)}{(2 Z_j)^{1/2}} \right) \right\} + n\frac{1 + 2 \hat{\varphi}_n}{4 \hat{\varphi}_n}}
#' with \eqn{\hat{\varphi}_n = \frac{\hat{\lambda}_n}{\hat{\mu}_n}}, where \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution. 
#' Furthermore \eqn{Z_{jk} = (Y_j + Y_k)} and \eqn{Z_j = Y_j}, where \eqn{Y_i = \frac{X_i}{\hat{\mu}_n}} for \eqn{(X_i)_{i = 1,...,n}}, a sequence of  independent observations of a nonnegative random variable \eqn{X}.
#' To ensure numerical stability of the implementation the exponentially scaled complementary error function \eqn{\text{erfce}(x)} is used: \eqn{\text{erfce}(x) = \exp{(x^2)}\text{erfc}(x)}, with \eqn{\text{erfc}(x) = 2\int_x^\infty \exp{(-t^2)}dt/\pi}.
#' The null hypothesis is rejected for large values of the test statistic \eqn{HK_{n,a}^{(2)}}.
#' @references
#' Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
#'
#' @examples
#' HK2(rmutil::rinvgauss(20,2,1))
#' 
#' @export
HK2 <- function(data){
  n       = length(data)
  muH     = mean(data)
  lambdaH = 1/mean(1/data-1/muH)
  phiH    = lambdaH/muH
  Y       = data/muH
  Tmat1   = matrix(0,n,n)
  for (j in 1:n){
    for (k in 1:n){
      Zjk        = Y[j]+Y[k]
      Tmat1[j,k] = 1/Zjk
    }
  }
  Tmat2 = rep(0,n)
  for (j in 1:n){
    Zj       = Y[j]
    T1       = sqrt(pi*phiH/(2*Zj))
    T2       = erfce(sqrt(phiH)*(Zj+1)/sqrt(2*Zj))
    Tmat2[j] = (1-T1*T2)/Zj
  }
  TestStat = sum(Tmat1)/n - 2*sum(Tmat2) + n*(1+2*phiH)/(4*phiH)
  return(TestStat)
}



#' The Kolmogorov-Smirnov test statistic
#'
#' @description
#' This function computes the test statistic of the goodness-of-fit test for the inverse Gaussian family in the spirit of Kolmogorov and Smirnov. Note that this tests the composite hypothesis of fit to the family of inverse Gaussian distributions.
#'
#' @param data  a vector of positive numbers.
#'
#' @return value of the test statistic.
#'
#' @details
#'  Let \eqn{X_{(j)}} denote the \eqn{j}th order statistic of \eqn{X_1, \ldots, X_n}, a sequence of  independent observations of a positive random variable \eqn{X}. Furthermore, let \eqn{\hat{F}(x) = F(x; \hat{\mu}_n, \hat{\lambda}_n)}, where \eqn{F} is the distribution function of the inverse Gaussian distribution. 
#'  Note that  \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution.
#'  The null hypothesis is rejected for large values of the test statistic: 
#'  \deqn{KS = \max(D^+, D^-),}
#'  where \deqn{D^+ = \max_{j=1,\ldots,n} \left( \frac{j}{n} - \hat{F}(X_{(j)}) \right)}
#'  and \deqn{D^- = \max_{j=1,\ldots,n} \left( \hat{F}(X_{(j)}) - \frac{j-1}{n} \right).}
#'  
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2022) "On Testing the Adequacy of the Inverse Gaussian Distribution". \href{https://www.mdpi.com/2227-7390/10/3/350}{LINK}
#' 
#' @examples
#' KS(rmutil::rinvgauss(20,2,1))
#' 
#' @export
KS <- function(data){
  n        = length(data)
  muH      = mean(data)
  lambdaH  = 1/mean(1/data-1/muH)
  Fhat     = rmutil::pinvgauss(sort(data),muH,1/lambdaH)
  Dplus    = max((1:n)/n-Fhat)
  Dminus   = max(Fhat-(0:(n-1))/n)
  TestStat = max(Dplus,Dminus)
  return(TestStat)
}



#' The Cramer-von Mises test statistic
#'
#' @description
#' This function computes value of the test statistic of the goodness-of-fit test for the inverse Gaussian family in the spirit of Cramer and von Mises. Note that this tests the composite hypothesis of fit to the family of inverse Gaussian distributions.
#'
#' @param data  a vector of positive numbers.
#'
#' @return value of the test statistic.
#' 
#' @details
#' Let \eqn{X_{(j)}} denote the \eqn{j}th order statistic of \eqn{X_1, \ldots, X_n}, a sequence of  independent observations of a positive random variable \eqn{X}. Furthermore, let \eqn{\hat{F}(x) = F(x; \hat{\mu}_n, \hat{\lambda}_n)}, where \eqn{F} is the distribution function of the inverse Gaussian distribution. 
#' Note that  \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution.
#' The null hypothesis is rejected for large values of the test statistic: 
#' \deqn{CM = \frac{1}{12n} + \sum_{j=1}^{n} \left( \hat{F}(X_{(j)}) - \frac{2j-1}{2n} \right)^2.}
#'  
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2022) "On Testing the Adequacy of the Inverse Gaussian Distribution". \href{https://www.mdpi.com/2227-7390/10/3/350}{LINK}
#' 
#' @examples
#' CM(rmutil::rinvgauss(20,2,1))
#' 
#' @export
CM <- function(data){
  n        = length(data)
  muH      = mean(data)
  lambdaH  = 1/mean(1/data-1/muH)
  Fhat     = rmutil::pinvgauss(sort(data),muH,1/lambdaH)
  TestStat = 1/(12*n)+sum((Fhat-(2*(1:n)-1)/(2*n))^2)
  return(TestStat)
}



#' The Anderson-Darling test statistic
#'
#' @description
#' This function computes the test statistic of the goodness-of-fit test for the inverse Gaussian family in the spirit of Anderson and Darling.
#'
#' @param data  a vector of positive numbers.
#'
#' @return value of the test statistic.
#'
#' @details
#' Let \eqn{X_{(j)}} denote the \eqn{j}th order statistic of \eqn{X_1, \ldots, X_n}, a sequence of  independent observations of a positive random variable \eqn{X}. Furthermore, let \eqn{\hat{F}(x) = F(x; \hat{\mu}_n, \hat{\lambda}_n)}, where \eqn{F} is the distribution function of the inverse Gaussian distribution. 
#' Note that  \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution.
#' The null hypothesis is rejected for large values of the test statistic: 
#' \deqn{AD = -n - \frac{1}{n} \sum_{j=1}^{n} \left[ (2j-1) \log \hat{F}(X_{(j)}) + (2(n-j) + 1) \log \left( 1 - \hat{F}(X_{(j)}) \right) \right].}
#' 
#'
#' @references
#' Allison, J.S., Betsch, S., Ebner, B., Visagie, I.J.H. (2022) "On Testing the Adequacy of the Inverse Gaussian Distribution". \href{https://www.mdpi.com/2227-7390/10/3/350}{LINK}
#' 
#' @examples
#' AD(rmutil::rinvgauss(20,2,1))
#'
#' @export
AD <- function(data){
  n        = length(data)
  muH      = mean(data)
  lambdaH  = 1/mean(1/data-1/muH)
  Fhat     = rmutil::pinvgauss(sort(data),muH,1/lambdaH)
  T1       = (2*(1:n)-1)*log(Fhat)
  T2       = (2*(n-(1:n))+1)*log(1-Fhat)
  sum1     = sum(T1+T2)
  TestStat = -n-sum1/n
  return(TestStat)
}


#' The Baringhaus-Gaigall test statistic
#'
#' @description
#' This function computes the test statistic of the goodness-of-fit test for the inverse Gaussian family due to Baringhaus and Gaigall (2015).
#'
#' @param data  a vector of positive numbers.
#'
#' @return value of the test statistic.
#'
#' @details 
#' The test statistic of the Baringhaus-Gaigall test is defined as: 
#' \deqn{BG_{n} = \frac{n}{(n(n-1))^5} \sum_{\mu, \nu = 1, \mu \neq \nu}^{n} \left( N_1(\mu, \nu)N_4(\mu, \nu) - N_2(\mu, \nu)N_3(\mu, \nu) \right)^2,}
#' where
#' \deqn{N_1(\mu, \nu) = \sum_{i,j = 1, i \neq j}^{n} \mathbf{1} \left\{ \tilde{Y}_{i,j} \leq \tilde{Y}_{\mu, \nu}, \tilde{Z}_{i,j} \leq \tilde{Z}_{\mu, \nu} \right\},}
#' \deqn{N_2(\mu, \nu) = \sum_{i,j = 1, i \neq j}^{n} \mathbf{1} \left\{ \tilde{Y}_{i,j} \leq \tilde{Y}_{\mu, \nu}, \tilde{Z}_{i,j} > \tilde{Z}_{\mu, \nu} \right\},}
#' \deqn{N_3(\mu, \nu) = \sum_{i,j = 1, i \neq j}^{n} \mathbf{1} \left\{ \tilde{Y}_{i,j} > \tilde{Y}_{\mu, \nu}, \tilde{Z}_{i,j} \leq \tilde{Z}_{\mu, \nu} \right\},}
#' \deqn{N_4(\mu, \nu) = \sum_{i,j = 1, i \neq j}^{n} \mathbf{1} \left\{ \tilde{Y}_{i,j} > \tilde{Y}_{\mu, \nu}, \tilde{Z}_{i,j} > \tilde{Z}_{\mu, \nu} \right\},}
#' with \eqn{\mathbf{1}} being the indicator function.
#' Let \eqn{f(X_i,X_j) = (X_i + X_j)/2} and \eqn{g(X_i,X_j) = (X_i^{-1} + X_j^{-1})/2 - f(X_i,X_j)^{-1}}, with \eqn{X_1,...,X_n} positive, independent and identically distributed random variables with finite moments \eqn{\mathbb{E}[X_1^2]} and \eqn{\mathbb{E}[X_1^{-1}]}. 
#' Then \eqn{(\tilde{Y}_{i,j}, \tilde{Z}_{i,j}) = (f(X_i,X_j), g(X_i,X_j)), 1 \leq i,j \leq n, i \neq j}. Note that \eqn{\tilde{Y}_{i,j}} and \eqn{\tilde{Z}_{i,j}} are independent if, and only if \eqn{X_1,...,X_n} are realized from an inverse Gaussian distribution.
#' 
#' @references
#' Baringhaus, L.  Gaigall, D. (2015). "On an independence test approach to the goodness-of-fit problem", Journal of Multivariate Analysis, 140, 193-208. \doi{https://doi.org/10.1016/j.jmva.2015.05.013}
#'
#' @examples
#' BG(rmutil::rinvgauss(20,2,1))
#' 
#' @export
BG <- function(data){
  n  = length(data)
  Xi = matrix(rep(data,n),n,n)
  Xj = t(matrix(rep(data,n),n,n))
  Y  = (Xi+Xj)/2
  Z  = (1/Xi+1/Xj)/2 - 1/Y
  Y  = matrix(Y[lower.tri(Y,diag=F)|upper.tri(Y,diag=F)],n-1,n)
  Z  = matrix(Z[lower.tri(Z,diag=F)|upper.tri(Z,diag=F)],n-1,n)

  N1 = matrix(0,n-1,n)
  N2 = matrix(0,n-1,n)
  N3 = matrix(0,n-1,n)
  N4 = matrix(0,n-1,n)

  for (u in 1:(n-1)){
    for (v in 1:n){
      N1[u,v] = sum((Y<=Y[u,v])*(Z<=Z[u,v]))
      N2[u,v] = sum((Y<=Y[u,v])*(Z>Z[u,v]))
      N3[u,v] = sum((Y>Y[u,v])*(Z<=Z[u,v]))
      N4[u,v] = sum((Y>Y[u,v])*(Z>Z[u,v]))
    }
  }

  TestStat = n/(n*(n-1))^5*sum((N1*N4-N2*N3)^2)
  return(TestStat)
}
