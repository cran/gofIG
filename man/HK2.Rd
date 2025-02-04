% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{HK2}
\alias{HK2}
\title{The second Henze-Klar test statistic}
\usage{
HK2(data)
}
\arguments{
\item{data}{a vector of positive numbers.}
}
\value{
value of the test statistic.
}
\description{
This function computes the test statistic of  the second goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
}
\details{
The representation of the second Henze-Klar test statistic used for computation \eqn{(a = 0)} is given by:
\deqn{HK_{n,0}^{(2)} = \frac{1}{n} \sum_{j,k=1}^{n} Z_{jk}^{-1} - 2 \sum_{j=1}^{n} Z_j^{-1} \left\{ 1 - \sqrt{\frac{\pi \hat{\varphi}_n}{2 Z_j}} \, \mathrm{erfce} \left( \frac{\hat{\varphi}_n^{1/2} (Z_j + 1)}{(2 Z_j)^{1/2}} \right) \right\} + n\frac{1 + 2 \hat{\varphi}_n}{4 \hat{\varphi}_n}}
with \eqn{\hat{\varphi}_n = \frac{\hat{\lambda}_n}{\hat{\mu}_n}}, where \eqn{\hat{\mu}_n,\hat{\lambda}_n} are the maximum likelihood estimators for \eqn{\mu} and \eqn{\lambda}, respectively, the parameters of the inverse Gaussian distribution. 
Furthermore \eqn{Z_{jk} = (Y_j + Y_k)} and \eqn{Z_j = Y_j}, where \eqn{Y_i = \frac{X_i}{\hat{\mu}_n}} for \eqn{(X_i)_{i = 1,...,n}}, a sequence of  independent observations of a nonnegative random variable \eqn{X}.
To ensure numerical stability of the implementation the exponentially scaled complementary error function \eqn{\text{erfce}(x)} is used: \eqn{\text{erfce}(x) = \exp{(x^2)}\text{erfc}(x)}, with \eqn{\text{erfc}(x) = 2\int_x^\infty \exp{(-t^2)}dt/\pi}.
The null hypothesis is rejected for large values of the test statistic \eqn{HK_{n,a}^{(2)}}.
}
\examples{
HK2(rmutil::rinvgauss(20,2,1))

}
\references{
Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
}
