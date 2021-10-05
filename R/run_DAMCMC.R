
#' Run a DA-MCMC to fit the stochastic SIR model to discretely observed incidence counts with the PD-SIR algorithm
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param Y observed data
#' @param theta_0 initial value for the parameters
#' @param print_i logical; whether to print the iteration
#' @param save_SS logical; whether to save the sufficient statistics generated each iteration
#' @param par_prior parameters of the prior distributions
#'
#' @return list with the draws for the parameters, the log likelihood, acceptance rate and size of initial susceptible population
#' @export
#'
run_DAMCMC <- function(
  Y, N = 1e4,
  rho = 1, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 1, theta_0,
  print_i = FALSE,  save_SS = FALSE,
  par_prior = list(
    a_beta = 0.01, b_beta = 1,
    a_gamma = 1, b_gamma = 1,
    a_R0 = 2, b_R0 = 2,
    a_lambda = 1, b_lambda = 1
  )
  ) {

  run_time <- system.time({

  # Setup
  SS_save <- x_save <- theta_save <- vector(mode = "list", length = N / thin)
  f_save  <- numeric(N / thin)
  accept  <- numeric(N)

  theta_0 <- complete_theta(theta_0, iota_dist, Y[["S0"]])
  theta_save[[1]] <- theta <- theta_0

  SS_current <- list(compatible = FALSE)
  while(! SS_current[["compatible"]]) {
    x_current  <- rprop_x(theta, Y, gener, b, iota_dist, approx)
    SS_current <- suff_stat(x_current, Y, gener, b)
  }

  f_save[1] <- f_target <- f_log(theta, SS_current, gener, b, iota_dist)

  #
  # MH
  for(i in 2 : N) {

    # Update theta (Gibbs)
    theta <- gibbs_theta(SS_current, iota_dist, par_prior, theta, param, Y)


    # Propose latent space
    x_new    <- rprop_x(theta, Y, gener, b, iota_dist, approx, x_current, rho)
    SS_new   <- suff_stat(x_new, Y, gener, b)
    i_update <- x_new[["i_update"]]

    # If x_new incompatible with observed data Y
    if(!SS_new[["compatible"]]) { # ... keep previous iteration.
      if(print_i)  print(paste(i, " - incompatible"))
      if(i %% thin == 0) {
        j <- i / thin
        theta_save[[j]] <- theta
        if(save_SS)  SS_save[[j]] <- SS_current
        f_save    [[j]] <- f_log(theta, SS_current, gener, b, iota_dist)
      }

      next
    }

    # Target density
    f_target_current <- f_log(theta, SS_current, gener, b, iota_dist)
    f_target_new     <- f_log(theta, SS_new    , gener, b, iota_dist)

    # Proposal density
    f_prop_current <- dprop_x(theta, Y, x_current, i_update, gener, b, iota_dist, approx)
    f_prop_new     <- dprop_x(theta, Y, x_new    , i_update, gener, b, iota_dist, approx)

    # MH ratio
    R_log     <- f_target_new - f_target_current - f_prop_new + f_prop_current
    R         <- min(1, exp(R_log))
    accept[i] <- stats::runif(1) < R

    # Accept/reject new draws
    if(accept[i]) {
      if(print_i)  print(i)
      x_current        <- x_new
      f_target_current <- f_target_new
      SS_current       <- SS_new
    } # end-if

    # Save thinned draws
    if(i %% thin == 0) {
      j <- i / thin
      theta_save[[j]] <- theta
      if(save_SS)  SS_save[[j]] <- SS_current
      f_save    [[j]] <- f_target_current
    }

  } # end-for

  # Output
  out <- list(
    theta = theta_save, loglik = f_save, rate_accept = mean(accept), S0 = Y[["S0"]]
  )
  if(save_SS) out[["SS"]] <- SS_save


  }) # end-system.time
  run_time <- as.numeric(run_time)[3]
  out[["run_time"]] <-  run_time

  return(out)

}
