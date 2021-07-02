

#' Gibbs sampler for theta
#'
#' allows two parameterizations of the stochastic SIR model
#'
#' @param SS sufficient statistics
#' @param par_prior parameters of the prior distribution
#' @param theta parameters of the SIR process
#' @param parameterization c("bg", "bR"); parameterize the models in terms of (beta, gamma) or (beta, R0)
#' @param Y observed data
#'
#' @return new value for theta
#' @export
#'
gibbs_theta <- function(SS, par_prior, theta, parameterization, Y){

  if(parameterization == "bg") { # theta = (beta, gamma)
    gibbs_beta_gamma(SS, par_prior, Y)
  } else if (parameterization == "bR") { # theta = (beta, R0)
    gibbs_beta_R0(SS, par_prior, theta, Y)
  }

}
