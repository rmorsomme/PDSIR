
# Proposal log density of the infection amd removal times tau_T and tau_J

#' Log proposal density of observed removals
#'
#' @param theta parameters of SIR model
#' @param tau_T infection times
#' @param tau_J removal times
#' @param tau_J_obs index of observed removal times
#' @param iota_distribution c("exponential", "weibull"); distribution of infectious time
#' @param t_end end of observation window
#'
#' @return log density
#' @export
#'
contribution_observed_removal <- function(
  theta, tau_T, tau_J, tau_J_obs, iota_distribution, t_end
) {

  # Setup
  gamma  <- theta[["gamma"]]
  lambda <- theta[["lambda"]]
  nu     <- theta[["nu"]]

  # loglik
  iota_obs <- tau_J[tau_J_obs] - tau_T[tau_J_obs]
  loglik <- if(iota_distribution == "exponential") {
    dexp_trunc_log    (iota_obs, gamma             , 0, t_end - tau_T[tau_J_obs])
  } else if(iota_distribution == "weibull") {
    dweibull2_trunc_log(iota_obs, nu, rate = lambda, 0, t_end - tau_T[tau_J_obs])
  }

  # Output
  return(loglik)

}
