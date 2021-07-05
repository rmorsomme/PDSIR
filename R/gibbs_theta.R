

#' Gibbs sampler for theta
#'
#' allows two parameterizations of the stochastic SIR model
#'
#' @param SS sufficient statistics
#' @param par_prior parameters of the prior distribution
#' @param theta parameters of the SIR process
#' @param param c("bg", "bR"); parameterize the models in terms of (beta, gamma) or (beta, R0)
#' @param Y observed data
#'
#' @return new value for theta
#' @export
#'
gibbs_theta <- function(SS, par_prior, theta, param, Y){

  # Setup
  n_T         <- SS[["n_T"        ]]
  n_J         <- SS[["n_J"        ]]
  integral_SI <- SS[["integral_SI"]]
  integral_I  <- SS[["integral_I" ]]

  a_beta      <- par_prior[["a_beta"]]
  b_beta      <- par_prior[["b_beta"]]
  a_gamma     <- par_prior[["a_gamma"]]
  b_gamma     <- par_prior[["b_gamma"]]
  a_R0        <- par_prior[["a_R0"  ]]
  b_R0        <- par_prior[["b_R0"  ]]

  beta_curr <- theta[["beta"]]
  R0_curr   <- theta[["R0"  ]]

  S0        <- Y[["S0"]]

  # Full conditionals
  if(param == "bg") { # theta = (beta, gamma)
    beta_new  <- stats::rgamma(1, a_beta  + n_T, b_beta  + integral_SI)
    gamma_new <- stats::rgamma(1, a_gamma + n_J, b_gamma + integral_I )
    R0_new    <- S0 * beta_new / gamma_new
  } else if (param == "bR") { # theta = (beta, R0)
    beta_new  <-     stats::rgamma(1, a_beta + n_J + n_T, b_beta + integral_SI + integral_I * S0 / R0_curr)
    R0_new    <- 1 / stats::rgamma(1, a_R0   + n_J      , b_R0   + beta_curr * S0 * integral_I            )
    gamma_new <- S0 * beta_new / R0_new
  }

  theta_new <- list(beta = beta_new, gamma = gamma_new, R0 = R0_new)
  return(theta_new)

}
