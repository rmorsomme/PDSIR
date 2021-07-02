
# allows non-Markovian dynamics through Weibull-distribution infectious periods (Streftaris and Gibson, 2004)
#
# Severo, N. C. (1969). Generalizations of some stochastic epidemic models. Mathematical Biosciences, 4(3-4), 395-402.
# Streftaris, G., & Gibson, G. J. (2004). Bayesian inference for stochastic epidemics in closed populations. Statistical Modelling, 4(1), 63-75.


#' Log likelihood of the stochastic SIR model
#'
#' @param theta parameters of the stochastic SIR process
#' @param SS sufficient statistics of some latent data
#' @param generalized logical; whether to use the generalized stochastic SIR process (Severo, 1969)
#' @param b parameter of the generalized SIR
#' @param iota_distribution c("exponential", "weibull"); distribution of the infection period
#'
#' @return log likelihood
#' @export
#'
f_log <- function(theta, SS, generalized, b, iota_distribution = "exponential") {

  # Setup
  beta   <- theta[["beta" ]]
  gamma  <- theta[["gamma"]]
  rate   <- theta[["rate" ]]
  shape  <- theta[["shape"]]

  n_T             <- SS[["n_T"]]
  I_tau_T         <- SS[["I_tau_T"]]
  S_tau_T         <- SS[["S_tau_T"]]
  integral_SI     <- SS[["integral_SI"]]
  iota_recovered  <- SS[["iota_recovered"]]
  iota_infectious <- SS[["iota_infectious"]]

  # # Loglik
  # loglik <- # agent-based likelihood
  #   n_T   * log(beta) + sum(log(I_tau_T)) +
  #   n_J   * log(gamma) -
  #   beta  * integral_SI -
  #   gamma * integral_I

  # contribution of infections
  loglik_infec <- if(generalized){
    n_T * log(beta) + sum(log(I_tau_T) - b * log(S_tau_T)) - beta  * integral_SI
  } else {
    n_T * log(beta) + sum(log(I_tau_T)                   ) - beta  * integral_SI
  }

  # contribution of removals
  loglik_remov <-
    if(iota_distribution == "exponential") {
      sum(stats::dexp(iota_recovered , gamma      , log   = TRUE                    )) + # removals before t_end (observed)
      sum(stats::pexp(iota_infectious, gamma      , log.p = TRUE, lower.tail = FALSE))   # removals after t_end
    } else if(iota_distribution == "weibull") {
      sum(dweibull2  (iota_recovered , shape, rate, log   = TRUE                    )) +
      sum(pweibull2  (iota_infectious, shape, rate, log.p = TRUE, lower.tail = FALSE))
    }

  loglik <- loglik_infec + loglik_remov

  return(loglik)

}
