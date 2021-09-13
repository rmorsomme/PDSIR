#' Random generator for the truncated weibull distribution
#'
#' shape and rate parameterization.
#'
#' @param n number of random values to generate
#' @param shape shape parameter
#' @param lambda rate parameter
#'
#' @family weibull2
#'
#' @return vector of n random variables following a truncated weibull distribution
#' @export
#'
rweibull2 <- function(n, shape, lambda) {

  #scale <- rate_to_scale(lambda, shape)
  #X     <- stats::rweibull(n, shape = shape, scale = scale)

  U <- stats::runif(n)
  X <- (- log(U) / lambda) ^ (1 / shape)

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
pweibull2 <- function(x, shape, lambda, log.p = FALSE, lower.tail = TRUE) {

  prob <- if(log.p) {
    if(!lower.tail) {
      - lambda * x^shape
    } else {
      log(1 - exp( - lambda * x^shape))
    }
  } else {
    if(!lower.tail) {
          exp( - lambda * x^shape)
    } else {
      1 - exp( - lambda * x^shape)
    }
  }

#  scale <- rate_to_scale(lambda, shape)
#  prob <- stats::pweibull(x, shape, scale, lower.tail, log.p)

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
dweibull2 <- function(x, shape, lambda, log = FALSE) {

  dens <- if(log) {
    log(lambda) + log(shape) + (shape-1) * log(x) - lambda * x^shape
  } else {
    lambda * shape * x^(shape-1) * exp( - lambda * x^shape)
  }


#  scale <- rate_to_scale(lambda, shape)
#  dens  <- stats::dweibull(x, shape, scale, log)

  return(dens)

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
rweibull2_trunc <- function(n, shape, lambda, lower, upper) {

  shift <-         exp(- lambda * lower ^ shape)
  const <- shift - exp(- lambda * upper ^ shape)
  #const <- pweibull2(upper, shape, lambda) - pweibull2(lower, shape, lambda) # sanity check
  #shift <- pweibull2(lower, shape, lambda, lower.tail = FALSE) # sanity check

  U <- stats::runif(n)
  X <- (- log(shift - const * U) / lambda) ^ (1 / shape)

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
dweibull2_trunc_log <- function(x, shape, lambda, lower, upper) {

  stopifnot(all(lower <= x), all(x <= upper))

  const   <- exp(- lambda * lower ^ shape) - exp(- lambda * upper ^ shape)
  logdens <- dweibull2(x, shape, lambda, log = TRUE) - log(const)



#  const  <- pweibull2(upper, shape, lambda) - pweibull2(lower, shape, lambda)
#  loglik <- dweibull2(x    , shape, lambda, log = TRUE) - log(const)

  return(logdens)

}
