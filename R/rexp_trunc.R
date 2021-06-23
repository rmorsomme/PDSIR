#' Random generator for the truncated exponential distribution
#'
#'
#' Generates random values from a truncated exponential distribution
#'
#' @param n integer
#' @param lambda positive double
#' @param lower positive double
#' @param upper positive double
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
