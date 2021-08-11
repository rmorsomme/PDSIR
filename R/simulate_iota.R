#' Simulate infection periods
#'
#' @inheritParams simulate_SEM
#' @param n number of infection periods to simulate
#' @param gamma parameter of the exponential distribution
#' @param shape parameter of the weibull distribution
#' @param rate parameter of the weibull distribution
#'
#' @return vector of random infectious durations
#' @export
simulate_iota <- function(n, iota_dist, gamma, shape, rate) {

  iota <- if(iota_dist == "exponential") {
    stats::rexp(n, gamma)
  } else if(iota_dist == "weibull") {
    rweibull2(n, shape, rate)
  }

  return(iota)

}

