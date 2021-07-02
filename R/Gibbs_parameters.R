
#' Gibbs sampler for beta and gamma
#'
#' @param SS sufficient statistics of some latent data
#' @param par_prior parameters of the prior distribution
#' @param Y observed data
#'
#' @return draws from the joint full conditional for beta and gamma
#' @export
#'
gibbs_beta_gamma <- function(SS, par_prior, Y) {

  # Setup
  n_T         <- SS[["n_T"]]
  n_J         <- SS[["n_J"]]
  integral_SI <- SS[["integral_SI"]]
  integral_I  <- SS[["integral_I"]]

  a_beta      <- par_prior[["a_beta" ]]
  b_beta      <- par_prior[["b_beta" ]]
  a_gamma     <- par_prior[["a_gamma"]]
  b_gamma     <- par_prior[["b_gamma"]]

  S0          <- Y[["S0"]]


  # Full conditional (independent)
  beta_new  <- stats::rgamma(1, a_beta  + n_T, b_beta  + integral_SI)
  gamma_new <- stats::rgamma(1, a_gamma + n_J, b_gamma + integral_I )

  theta_new <- list(beta = beta_new, gamma = gamma_new, R0 = S0 * beta_new / gamma_new)
  return(theta_new)

}
