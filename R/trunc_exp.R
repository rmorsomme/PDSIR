
#' Random generator for the truncated exponential distribution
#'
#' Generates random values from a truncated exponential distribution
#'
#' @inheritParams propose_tau_T
#'
#' @param n number of random values to generates
#' @param lambda rate parameter
#'
#' @family truncated exponential
#'
#' @return positive double
#' @export
#'
#' @examples
#' hist(
#' rexp_trunc(1e4, lambda = 3, lower = 1, upper = 2),
#' main = "Truncated exponential distribution", xlab = NULL
#' )
#'
rexp_trunc <- function(n, lambda, lower, upper){ # Random Generator

  const <- exp(- lambda * lower) - exp(- lambda * upper) # normalizing constant
  # const == pexp(u, lambda) - pexp(l, lambda) # sanity check
  shift <- exp(- lambda * lower)
  U     <- stats::runif(n)
  X     <- - log(shift - U * const) / lambda
  return(X)

}


#' Log density of the truncated exponential distribution
#'
#' @inheritParams rexp_trunc
#'
#' @param x vector of values
#'
#' @family truncated exponential
#'
#' @return log density
#' @export
#'
dexp_trunc_log <- function(x, lambda, lower, upper){

  stopifnot(all(lower <= x), all(x <= upper))

  const  <- exp(- lambda * lower) - exp(- lambda * upper)
  loglik <- log(lambda) - lambda * x - log(const)

  return(loglik)

}
