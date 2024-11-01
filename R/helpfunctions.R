# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~ Our first test: subroutines ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# First subroutine required
h1 <- function(X1,X2,a){
  if(X1<=X2){
    return((a^-3)*(2-exp(-a*X1)*((a*X1+1)^2+1))+(a^-2)*X1*(exp(-a*X1)*(a*X1+1)-exp(-a*X2)*(a*X2+1))+(1/a)*X1*X2*exp(-a*X2))
  } else {
    return((a^-3)*(2-exp(-a*X2)*((a*X2+1)^2+1))+(a^-2)*X2*(exp(-a*X2)*(a*X2+1)-exp(-a*X1)*(a*X1+1))+(1/a)*X1*X2*exp(-a*X1))
  }
}

# Second subroutine required
h2 <- function (X1,X2,a){
  if(X1<=X2){
    return(X1*exp(-a*X2)/a)
  } else {
    return((a^-2)*(exp(-a*X2)*(a*X2+1)-exp(-a*X1)*(a*X1+1))+X1*exp(-a*X1)/a)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ Our second test: subroutines ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h1t <- function(s,t,a){
  if(s<=t){
    T1 = sqrt(pi/(a^3))/4
    T2 = s/(2*a)*exp(-a*s^2)
    T3 = a/2*sqrt(pi/(a^5))*stats::pnorm(-sqrt(2*a)*s)
    T4 = s/(2*a)*(exp(-a*s^2)-exp(-a*t^2))
    T5 = s*t*sqrt(pi/a)*stats::pnorm(-sqrt(2*a)*t)
    return(T1-T2-T3+T4+T5)
  } else {
    T1 = sqrt(pi/(a^3))/4
    T2 = t/(2*a)*exp(-a*t^2)
    T3 = a/2*sqrt(pi/(a^5))*stats::pnorm(-sqrt(2*a)*t)
    T4 = t/(2*a)*(exp(-a*t^2)-exp(-a*s^2))
    T5 = s*t*sqrt(pi/a)*stats::pnorm(-sqrt(2*a)*s)
    return(T1-T2-T3+T4+T5)
  }
}

h2t <- function (s,t,a){
  if(s<=t){
    return(s*sqrt(pi/a)*stats::pnorm(-sqrt(2*a)*t))
  } else {
    T1 = (exp(-a*t^2)-exp(-a*s^2))/(2*a)
    T2 = s*sqrt(pi/a)*stats::pnorm(-sqrt(2*a)*s)
    return(T1+T2)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~ Henze Klar (denoted Tna in Henze & Klar) ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Subroutine
erfce<-function(x){
  exp(x^2)*pracma::erfc(x)
}


#' Print method for tests of the inverse Gaussian distribution
#' @description
#' Printing objects of class "gofIG".
#'
#' @param x object of class "gofIG".
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' A \code{gofIG} object is a named list of numbers and character string, supplemented with \code{test} (the name of the teststatistic). \code{test} is displayed as a title.
#' The remaining elements are given in an aligned "name = value" format.
#'
#' @return
#' the argument x, invisibly, as for all \link{print} methods.
#'
#' @examples
#' \donttest{print(test.ABEV1(rgamma(20,1)))}
#'
#' @export
print.gofIG <- function(x, ...){
  cat("\n \n")
  cat("------------------------------------------------------------------------- \n")
  cat("\n")
  cat("         One-sample test for inverse Gaussian distribution with the", x$Test, " teststatistic.\n"  )
  cat("\n")
  if(!is.null(x$parameter)){
  cat("tuning parameter            = ", x$parameter," \n")
  }
  cat("value of statistic          = ", x$T.value, " \n")
  cat("p value                     = ", x$p.value, " \n")
  if(!is.null(x$est.method)){
  cat("estimation method           = ", x$est.method, " \n")
  }
  cat("estimated parameters        = ", x$par.est, " \n")
  cat("number of bootstrap samples = ", x$boot.run, " \n \n")
  cat("\n")
  cat("\n")
  cat("------------------------------------------------------------------------- \n")
}
