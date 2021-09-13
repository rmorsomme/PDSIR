#' Log proposal density of observed removals
#'
#' @inheritParams dprop_x
#'
#' @param theta parameters of SIR model
#' @param tau_T infection times
#' @param tau_J removal times
#' @param t_end end of observation window
#' @param i_update_obs indices of observations to update
#'
#' @return log density
#' @export
#'
contribution_observed_removal <- function(
  theta, tau_T, tau_J, i_update_obs, iota_dist, t_end
) {

  # Setup
  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  # loglik
  iota_obs <- tau_J[i_update_obs] - tau_T[i_update_obs]
  loglik <- if(iota_dist == "exponential") {
    dexp_trunc_log    (iota_obs, gamma         , 0, t_end - tau_T[i_update_obs])
  } else if(iota_dist == "weibull") {
    dweibull2_trunc_log(iota_obs, shape, lambda, 0, t_end - tau_T[i_update_obs])
  }

  # Output
  return(loglik)

}
