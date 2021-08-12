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
#' @param param c("bg", "bR"); parameterize the models in terms of (beta, gamma) or (beta, R0)
#' @param approx c("poisson", "ldp"); whether to approximate the distribution of the infection times with a poisson process or a linear death process
#' @param iota_dist c("exponential", "weibull"); distribution of the infection period
#' @param gener logical; whether to use the generalized SIR of Severo (1972)
#' @param b parameter of the generalized SIR
#' @param par_prior parameters of the prior distribution
#' @param theta_0_factor factors by which the true value of the parameters is multiplied to initialize the Markov chain
#'
#' @return list containing the parameters, observed data, Markov chain and run time of the algorithm
#' @export
#'
experiment_1_proof_of_concept <- function(
  S0 = 1e3, I0 = 1e1, theta = list(R0 = 2.5, gamma = 1),
  t_end = 6, K = 10,
  N = 1e4, thin = 1, rho = 1,
  param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  theta_0_factor = c(1, 1), # for c(beta, gamma)
  par_prior = list(a_beta = 0.1, b_beta = 1, a_gamma = 1, b_gamma = 1, a_R0 = 2, b_R0 = 2e-3)
  ) {

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]
  beta  <- R0 / gamma / S0

  SIR   <- simulate_SEM(S0, I0, t_end, theta)
  Y     <- observed_data(SIR, K)

  # run a long chain
  theta_0 <- list(
    beta  = theta_0_factor[1]  * beta ,
    gamma = theta_0_factor[2] * gamma,
    R0    = S0 * theta_0_factor[1] * beta / (theta_0_factor[2] * gamma)
    )
  {
    t0 <- Sys.time()
    MC <- run_DAMCMC(
      Y, N,
      rho, param, approx, par_prior,
      iota_dist, gener,
      thin, theta_0 = theta_0
      )
    run_time <- Sys.time() - t0
  }

  return(list(theta = theta, Y = Y, MC = MC, run_time = run_time))

}



#' Analyze output of MCMC run
#'
#' @param x object returned by the function experiment_1_proof_of_concept
#' @param plot_id name file for the figures
#' @param burnin number of iterations to discard, the default value NULL will discard half of the draws.
#' @param path directory in which to save the figure
#'
#' @return list of summary statistics with and without burnin
#' @export
#'
experiment_1_output_analysis <- function(
  x, plot_id = NULL, path = NULL, burnin = 0
  ) {

  theta <- x    [["theta"]]
  Y     <- x    [["Y"    ]]
  MC    <- x    [["MC"   ]]
  gamma <- theta[["gamma"]]
  beta  <- theta[["beta" ]]

  summary_no_burn <- analyze_MCMC(
    MC, burnin = 0,
    #beta_true = beta, gamma_true = gamma,
    plot_id = paste0(plot_id, "_no_burn"), path = path
    )

  if(is.null(burnin))  burnin <- length(MC[["theta"]])/2
  summary_burn <- analyze_MCMC(
    MC, burnin = burnin,
    beta_true = beta, gamma_true = gamma,
    plot_id = paste0(plot_id, "_burn"), path = path
    )

  out <- list(summary_no_burn = summary_no_burn, summary_burn = summary_burn)
  return(out)

}




#' Experiment 2: Compare trajectories of SIR and PD-SIR process
#'
#' @inheritParams experiment_1_proof_of_concept
#' @inheritParams experiment_1_output_analysis
#'
#' @param Ks vector of numbers of time intervals to consider for the PD-SIR process
#'
#' @return figures comparing the trajectories of the PDSIR and SIR for different value K
#' @export
#'
experiment_2_PDSIR_trajectories <- function(
  S0 = 1e4, I0 = 1e1, theta = list(R0 = 4, gamma = 1),
  t_end = 6, Ks = c(3, 5, 10, 50, 1e3),
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  path
  ) {

  theta <- add_beta(theta, S0)

  # SIR
  SIR   <- simulate_SEM(S0, I0, t_end, theta)
  draw_trajectories(SIR, plot_id = "E2", path, t_end, type = "SIR")

  # PDSIR
  for(K in Ks) {

    Y        <- observed_data(SIR, K)
    PDSIR    <- rprop_x(theta, Y, gener, b, iota_dist, approx = "ldp")
    PDSIR_SI <- suff_stat(PDSIR, Y, gener, b, return_SI = TRUE)
    compare_trajectories(SIR, PDSIR_SI, paste0("E2_K", K), path, t_end)

  }

}
