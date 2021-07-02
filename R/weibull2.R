#' CDF of weibull distribution
#'
#' @param x vector
#' @param shape shape parameter
#' @param rate rate parameter
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x), otherwise, P(X > x)
#'
#' @return vector of tail probabilities
#' @export
#'
pweibull2 <- function(x, shape, rate, log.p = FALSE, lower.tail = TRUE) {

  scale <- rate_to_scale(rate, shape)
  stats::pweibull(x, shape, scale = scale, log.p = log.p, lower.tail = lower.tail)

}


#' Random generator for the truncated weibull distribution
#'
#' @param n number of random variables to generate.
#' @param shape shape parameter
#' @param rate rate parameter
#' @param lower lower bound
#' @param upper upper bound
#'
#' @return n-length vector of random variables following the weibull distribution
#' @export
#'
rweibull2_trunc <- function(n, shape, rate, lower, upper) {

  const <- pweibull2(upper, shape, rate) - pweibull2(lower, shape, rate)
  shift <- pweibull2(lower, shape, rate)

  U <- stats::runif(n)
  V <- U * const + shift
  X <- (- log(1 - V) / rate) ^ (1 / shape)

  return(X)

}


#' Density of weibull distribution
#'
#' @param x vector
#' @param shape shape parameter
#' @param rate rate parameter
#' @param log logical; whether to take the logarithm of the density
#'
#' @return vector of densities
#' @export
#'
dweibull2 <- function(x, shape, rate, log = FALSE) {

  scale <- rate_to_scale(rate, shape)
  stats::dweibull(x, shape, scale = scale, log = log)

}


#

#' Log density of the truncated Weibull distribution
#'
#' @param x vector
#' @param shape shape parameter
#' @param rate rate parameter
#' @param lower lower bound
#' @param upper uppaer bound
#'
#' @return vector of densities
#' @export
#'
dweibull2_trunc_log <- function(x, shape, rate, lower, upper) {

  stopifnot(all(lower <= x), all(x <= upper))

  const  <- pweibull2(upper, shape, rate) - pweibull2(lower, shape, rate)
  loglik <- dweibull2(x    , shape, rate, log = TRUE) - log(const)

  return(loglik)

}
