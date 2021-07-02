
#' Generates removal times
#'
#' @param theta parameters of the SIR process
#' @param tau_T infection times
#' @param iota_distribution c("exponential", "weibull"); distribution of infection period
#' @param t_end end of observation period
#'
#' @return removal times
#' @export
#'
propose_tau_J <- function(theta, tau_T, iota_distribution, t_end) {

  # Setup
  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  nu     <- theta[["nu"    ]]

  # Compute p_i
  iota_obs <- if(iota_distribution == "exponential") {
    rexp_trunc    (length(tau_T), gamma             , 0, t_end - tau_T)
  } else      if(iota_distribution == "weibull") {
    rweibull2_trunc(length(tau_T), nu, rate = lambda, 0, t_end - tau_T)
  }

  tau_J <- tau_T + iota_obs

  # Output
  return(tau_J)

}
