

#
#' Log density of the truncated exponential distribution
#'
#' @param x vector
#' @param lambda rate of exponential distribution
#' @param lower lower bound
#' @param upper upper bound
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
