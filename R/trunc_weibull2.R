#' Random generator for the truncated weibull distribution
#'
#' shape and rate parameterization.
#'
#' @param n number of random values to generate
#' @param shape shape parameter
#' @param rate rate parameter
#'
#' @family weibull2
#'
#' @return vector of n random variables following a truncated weibull distribution
#' @export
#'
rweibull2 <- function(n, shape, rate) {


  scale <- rate_to_scale(rate, shape)
  X     <- stats::rweibull(n, shape = shape, scale = scale)

  return(X)

}

#' CDF of weibull distribution
#'
#' Lower and upper tail probabilities of the weibull distribution with the shape and rate parameterization.
#'
#' @inheritParams rweibull2
#' @param x vector of values
#' @param log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x), otherwise, P(X > x)
#'
#' @family weibull2
#'
#' @return vector of tail probabilities
#' @export
#'
pweibull2 <- function(x, shape, rate, log.p = FALSE, lower.tail = TRUE) {

  scale <- rate_to_scale(rate, shape)
  prob <- stats::pweibull(
    x, shape = shape, scale = scale, log.p = log.p, lower.tail = lower.tail
    )

  return(prob)

}

#' Density of weibull distribution
#'
#' @inheritParams pweibull2
#' @param log logical; whether to take the logarithm of the density
#'
#' @family weibull2
#'
#' @return vector of densities
#' @export
#'
dweibull2 <- function(x, shape, rate, log = FALSE) {

  scale <- rate_to_scale(rate, shape)
  prob <- stats::dweibull(x, shape, scale = scale, log = log)

  return(prob)

}





#' Random generator for the truncated weibull distribution
#'
#' @inheritParams propose_tau_T
#' @inheritParams pweibull2
#'
#' @family weibull2
#'
#' @return vector of n random variables following a truncated weibull distribution
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


#' Log density of the truncated Weibull distribution
#'
#' @inheritParams rweibull2_trunc
#' @param x vector of values
#'
#' @family weibull2
#'
#' @return vector of log densities
#' @export
#'
dweibull2_trunc_log <- function(x, shape, rate, lower, upper) {

  stopifnot(all(lower <= x), all(x <= upper))

  const  <- pweibull2(upper, shape, rate) - pweibull2(lower, shape, rate)
  loglik <- dweibull2(x    , shape, rate, log = TRUE) - log(const)

  return(loglik)

}
