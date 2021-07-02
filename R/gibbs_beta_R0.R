
#' Two-stage Gibbs sampler for beta and R0 for the stochastic SIR model
#'
#' @param SS sufficient statistics#'
#' @param par_prior parameters of the prior distribution
#' @param theta parameters of the SIR process
#' @param Y observed data
#'
#' @return new value for (beta, R0)
#' @export
#'
gibbs_beta_R0 <- function(SS, par_prior, theta, Y) {

  # Setup
  n_T         <- SS[["n_T"        ]]
  n_J         <- SS[["n_J"        ]]
  integral_SI <- SS[["integral_SI"]]
  integral_I  <- SS[["integral_I" ]]

  a_beta      <- par_prior[["a_beta"]]
  b_beta      <- par_prior[["b_beta"]]
  a_R0        <- par_prior[["a_R0"  ]]
  b_R0        <- par_prior[["b_R0"  ]]

  beta_curr <- theta[["beta"]]
  R0_curr   <- theta[["R0"  ]]

  S0        <- Y[["S0"]]

  # Full conditional distributions
  beta_new  <-     stats::rgamma(1, a_beta + n_J + n_T, b_beta + integral_SI + integral_I * S0 / R0_curr)
  R0_new    <- 1 / stats::rgamma(1, a_R0   + n_J      , b_R0   + beta_curr * S0 * integral_I            )
  gamma_new <- S0 * beta_new / R0_new

  theta_new <- list(beta = beta_new, gamma = gamma_new, R0 = R0_new)
  return(theta_new)

}
