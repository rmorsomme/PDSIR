#' Gibbs sampler for theta
#'
#' allows two parameterizations of the stochastic SIR model
#'
#' @inheritParams run_DAMCMC
#'
#' @param SS sufficient statistics
#' @param par_prior parameters of the prior distributions
#' @param theta parameters of the SIR process
#' @param param c("bg", "bR"); parameterization (beta, gamma) or (beta, R0)
#' @param Y observed data
#'
#' @return new value for theta
#' @export
#'
gibbs_theta <- function(SS, iota_dist, par_prior, theta, param, Y){

  # Setup
  n_T             <- SS[["n_T"            ]]
  n_J             <- SS[["n_J"            ]]
  integral_SI     <- SS[["integral_SI"    ]]
  integral_I      <- SS[["integral_I"     ]]
  iota_infectious <- SS[["iota_infectious"]]
  iota_removed    <- SS[["iota_removed"   ]]

  a_beta   <- par_prior[["a_beta"  ]]
  b_beta   <- par_prior[["b_beta"  ]]
  a_gamma  <- par_prior[["a_gamma" ]]
  b_gamma  <- par_prior[["b_gamma" ]]
  a_R0     <- par_prior[["a_R0"    ]]
  b_R0     <- par_prior[["b_R0"    ]]
  a_lambda <- par_prior[["a_lambda"]]
  b_lambda <- par_prior[["b_lambda"]]

  beta_curr   <- theta[["beta"  ]]
  R0_curr     <- theta[["R0"    ]]
  lambda_curr <- theta[["lambda"]]
  shape       <- theta[["shape" ]]

  S0          <- Y[["S0"]]

  # Full conditionals (Gibbs)

  if(iota_dist == "exponential") { # Markovian process - 2 parameterizations

    if(param == "bg") { # theta = (beta, gamma)

      beta_new   <- stats::rgamma(1, a_beta  + n_T, b_beta  + integral_SI)
      gamma_new  <- stats::rgamma(1, a_gamma + n_J, b_gamma + integral_I )
      R0_new     <- S0 * beta_new / gamma_new

    } else if (param == "bR") { # theta = (beta, R0)

      beta_new   <-     stats::rgamma(1, a_beta + n_J + n_T, b_beta + integral_SI + integral_I * S0 / R0_curr)
      R0_new     <- 1 / stats::rgamma(1, a_R0   + n_J      , b_R0   + beta_new * S0 * integral_I             )
      gamma_new  <- S0 * beta_new / R0_new

    } # end-if parameterization

    theta_new  <- list(R0 = R0_new, gamma = gamma_new, beta = beta_new)

  } else if(iota_dist == "weibull") { # non-Markovian process

    # Update beta
    beta_new  <- stats::rgamma(1, a_beta + n_T, b_beta + integral_SI)

    # Update lambda
    sum_iota   <- sum(iota_infectious ^ shape) + sum(iota_removed ^ shape)
    lambda_new <- stats::rgamma(1, a_lambda + n_J, b_lambda + sum_iota)

    # Complete theta
    gamma_new  <- 1 / (lambda_new ^ (- 1 / shape) * gamma(1 + 1 / shape))
    R0_new     <- S0 * beta_new / gamma_new

    theta_new  <- list(
      R0 = R0_new, lambda = lambda_new, shape = shape, gamma = gamma_new, beta = beta_new
      )

  } # end-if iota_dist

  return(theta_new)

}


