#' Simulate infection periods
#'
#' @inheritParams simulate_SEM
#' @param n number of infection periods to simulate
#' @param gamma parameter of the exponential distribution
#' @param shape parameter of the weibull distribution
#' @param lambda parameter of the weibull distribution
#'
#' @return vector of random infectious durations
#' @export
simulate_iota <- function(n, iota_dist, gamma, shape, lambda) {

  iota <- if(iota_dist == "exponential") {
   # stats::rexp(n, gamma)
    rexp2(n, gamma)
  } else if(iota_dist == "weibull") {
    rweibull2(n, shape, lambda)
  }

  return(iota)

}

