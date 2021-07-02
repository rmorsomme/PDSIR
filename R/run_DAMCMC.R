
#' Run a DA-MCMC to fit the stochastic SIR model to discretely observed incidence counts with the PD-SIR algorithm
#'
#' @param Y observed data
#' @param rho portion of the latent being updated each iteration
#' @param N number of iterations of the Markov chain
#' @param thin thinning argument for the Markov chain
#' @param theta_0 initial value for the parameters
#' @param par_prior parameters of the prior distribution
#' @param parameterization c("bg", "bR"); parameterize the models in terms of (beta, gamma) or (beta, R0)
#' @param generalized logical; whether to use the generalized SIR
#' @param b parameter of the generalized SIR
#' @param print_i logical; whether to print the iteration
#' @param save_SS logical; whether to save the sufficient statistics generated each iteration
#'
#' @return list with the draws for the parameters, the log likelihood, acceptance rate and size of initial susceptible population
#' @export
#'
run_DAMCMC <- function(
  Y, rho = 1/10, N = 1e4, thin = 1, theta_0,
  par_prior = list(a_beta = 0.1, b_beta = 1, a_gamma = 1, b_gamma = 1, a_R0 = 2, b_R0 = 2e-3),
  parameterization = "bR",
  generalized = FALSE, b = 1/2,
  print_i = FALSE, save_SS = FALSE
) {


  # Setup
  SS_save <- x_save <- theta_save <- vector(mode = "list", length = N / thin)
  f_save  <- numeric(N / thin)
  accept  <- numeric(N)

  theta <- theta_0
  SS_current <- list(compatible = FALSE)
  while(! SS_current[["compatible"]]) {
    x_current  <- rprop_x(theta, Y, generalized, b)
    SS_current <- sufficient_statistics(x_current, Y, generalized, b)
  }

  f_target   <- f_log(theta, SS_current, generalized, b)

  #
  # MH
  for(i in 2 : N) {

    # Propose theta (Gibbs)
    theta <- gibbs_theta(SS_current, par_prior, theta, parameterization, Y)


    # Propose latent space
    x_new    <- rprop_x(theta, Y, generalized, b, x_current, rho)
    SS_new   <- sufficient_statistics(x_new, Y, generalized, b)
    i_update <- x_new[["i_update"]]

    # If x_new incompatible with observed data Y
    if(!SS_new[["compatible"]]) { # ... keep previous iteration.
      if(print_i)  print(paste(i, " - incompatible"))
      if(i %% thin == 0) {
        j <- i / thin
        theta_save[[j]] <- theta
        SS_save   [[j]] <- SS_current
        f_save    [[j]] <- f_log(theta, SS_current, generalized, b)
      }

      next
    }

    # Target density
    f_target_current <- f_log(theta, SS_current, generalized, b)
    f_target_new     <- f_log(theta, SS_new    , generalized, b)

    # Proposal density
    f_prop_current <- dprop_x(
      theta, Y, x_current, i_update = i_update, generalized, b
    )
    f_prop_new     <- dprop_x(
      theta, Y, x_new    , i_update = i_update, generalized, b
    )

    # MH ratio
    R_log     <- f_target_new - f_target_current - f_prop_new + f_prop_current
    R         <- min(1, exp(R_log))
    accept[i] <- stats::runif(1) < R

    # Accept/reject new draws
    if(accept[i]) {
      if(print_i)  print(i)
      x_current          <- x_new
      f_target_current   <- f_target_new
      SS_current         <- SS_new
    } # end-if

    # Save thinned draws
    if(i %% thin == 0) {
      j <- i / thin
      theta_save[[j]] <- theta
      SS_save   [[j]] <- SS_current
      f_save    [[j]] <- f_target_current
    }

  } # end-for

  # Output
  out <- list(
    theta = theta_save, loglik = f_save, rate_accept = mean(accept), S0 = Y[["S0"]]
  )
  if(!save_SS) out[["SS"]] <- SS_save
  return(out)

}
