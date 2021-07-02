
#' Experiment 1: Proof of concept
#'
#' @param S0 size of initial susceptible population
#' @param I0 size initial infectious population
#' @param theta parameters of the SIR process
#' @param t_end end of observation period
#' @param K number of observation intervals
#' @param N number of iterations of the Markov chain
#' @param thin thinning parameter for the Markov chain
#' @param rho proportion of the latent data updated each iteration
#' @param theta_0_factor_beta factor by which the true value of beta is multiplied to determine the starting point of the Markov chain
#' @param theta_0_factor_gamma factor by which the true value of beta is multiplied to determine the starting point of the Markov chain
#'
#' @return list containing the parameters, observed data, Markov chain and run time of the algorithm
#' @export
#'
experiment_1_proof_of_concept <- function(
  S0 = 1e3, I0 = 1e1, theta = list(R0 = 2, gamma = 1),
  t_end = 6, K = 10,
  N = 1e5, thin = 1,
  rho = 1/5,
  theta_0_factor_beta = 0.1, theta_0_factor_gamma = 0.1
  ) {

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]
  beta  <- R0 / gamma / S0

  SIR  <- simulate_SIR(S0, I0, beta, gamma, t_end)
  Y    <- observed_data(SIR, K)

  # run a long chain
  t0 <- Sys.time()
  MC <- run_DAMCMC(
    Y, rho, N, thin, parameterization = "bR",
      theta_0 = list(
        beta  = theta_0_factor_beta  * beta ,
        gamma = theta_0_factor_gamma * gamma,
        R0 = S0 * theta_0_factor_beta * beta / (theta_0_factor_gamma * gamma)
        )
    )
    run_time <- Sys.time() - t0

  return(list(theta = theta, Y = Y, MC = MC, run_time = run_time))

}
