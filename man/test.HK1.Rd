% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IGtest.R
\name{test.HK1}
\alias{test.HK1}
\title{The first Henze-Klar goodness-of-fit test for the inverse Gaussian family}
\usage{
test.HK1(data, a = 0, B = 500)
}
\arguments{
\item{data}{a vector of positive numbers.}

\item{a}{positive tuning parameter.}

\item{B}{number of bootstrap iterations used to obtain p value.}
}
\value{
a list containing the value of the name of the test statistic, the used tuning parameter, the parameter estimation method, the value of the test statistic, the bootstrap p value, the values of the estimators, and the number of bootstrap iterations: \cr
\describe{
        \item{\code{$Test}}{the name of the used test statistic.}
        \item{\code{$parameter}}{the value of the tuning parameter.}
        \item{\code{$T.value}}{the value of the test statistic.}
        \item{\code{$p.value}}{the approximated p value.}
        \item{\code{$par.est}}{the estimated parameters.}
        \item{\code{$boot.run}}{number of bootstrap iterations.}
}
}
\description{
This function computes the goodness-of-fit test for the inverse Gaussian family due to Henze and Klar (2002).
}
\details{
The test statistics is a weighted integral over the squared modulus of some measure of deviation of the empirical distribution of given data from the family of inverse Gaussian laws, expressed by means of the empirical Laplace transform.
}
\examples{
test.HK1(rmutil::rinvgauss(20,2,1),B=100)

}
\references{
Henze, N. and Klar, B. (2002) "Goodness-of-fit tests for the inverse Gaussian distribution based on the empirical Laplace transform", Annals of the Institute of Statistical Mathematics, 54(2):425-444. \doi{https://doi.org/10.1023/A:1022442506681}
}
