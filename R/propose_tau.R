
#' Generates infection times
#'
#' @inheritParams rprop_x
#'
#' @param n number of infection times to generate
#' @param mu individual rate of infection
#' @param lower lower bound
#' @param upper upper bound
#'
#' @return infection times
#' @export
#'
propose_tau_T <- function(n, mu, lower, upper, approx) {

  stopifnot(approx %in% c("poisson", "ldp"))

  tau_T <- if(approx == "poisson") {
    stats::runif(n,     lower, upper)
  } else if(approx == "ldp") {
    rexp_trunc  (n, mu, lower, upper)
  }

  return(tau_T)

}



#' Generates removal times
#'
#' @inheritParams compute_pi
#'
#' @return removal times
#' @export
#'
propose_tau_J <- function(theta, tau_T, iota_dist, t_end) {

  # Setup
  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  # Compute p_i
  iota_obs <- if(iota_dist == "exponential") {
    rexp_trunc    (length(tau_T), gamma         , 0, t_end - tau_T)
  } else      if(iota_dist == "weibull") {
    rweibull2_trunc(length(tau_T), shape, lambda, 0, t_end - tau_T)
  }

  tau_J <- tau_T + iota_obs

  # Output
  return(tau_J)

}

