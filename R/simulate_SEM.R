#' Simulate trajectories of the Stochastic SIR and SEIR processes
#'
#' Also computes additional statistics such as MLEs and summary statistics
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param E0 intial exposed population size
#' @param type c("SIR", "SEIR"); type of process to simulate
#'
#' @return a list with all types of useful objects
#' @export
#'
simulate_SEM <- function(
  S0 = 1e3, E0 = 0, I0 = 1e1, t_end = 6,
  theta = list(R0 = 2, gamma = 1, shape = 2, rate = 1, epsilon = 1),
  gener = FALSE, b = 1/2,
  type = "SIR", # "SEIR"
  iota_dist = "exponential" # "weibull"
) {

  #
  # Population

  N0 <- S0 + I0


  #
  # Parameters

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]
  shape <- theta[["shape"]]
  rate  <- theta[["rate" ]]
  delta <- theta[["delta"]]

  iota_mean <- if(iota_dist == "exponential") {
    1 / gamma
  } else if(iota_dist == "weibull") {
    rate ^ (- 1 / shape) * gamma(1 + 1 / shape) # mean infectious time
  }

  beta  <- R0 / iota_mean / S0


  #
  # Initialization

  tau_T <- tau_F <- tau_J <- rep(Inf, N0)
  S <- E <- I <- t <- W <- X <- I_tau_t_true <- c() # S, E, I, R, time, waiting time, type of event, and number of infectious at infection times

  # Initialize tau's
  if(type == "SIR") {

    tau_T[1 : I0] <- 0 # TODO: relax assumption that individual initially infected at 0; important for non-Markovian process
    iotas <- simulate_iota(I0, iota_dist, gamma, shape, rate)
    tau_J[1 : I0] <- tau_T[1 : I0] + iotas

  } else if(type == "SEIR") {

    if(I0 > 0) {
      tau_T[1 : I0] <- -Inf
      tau_F[1 : I0] <- 0 # TODO: relax assumption
      iotas <- simulate_iota(I0, iota_dist, gamma, shape, rate)
      tau_J[1 : I0] <- tau_F[1 : I0] + iotas
    } # end-if(I0)
    if(E0 > 0) {
      tau_T[I0 + (1 : E0)] <- 0
      epsilons <- stats::rexp(E0, delta)
      tau_F[I0 + (1 : E0)] <- tau_T[I0 + (1 : E0)] + epsilons
    } # end-if(E0)
  }# end-if(SIR)

  # Initialize compartments, time, event type and event number
  S[1] <- S0
  E[1] <- E0
  I[1] <- I0
  t[1] <- 0
  X[1] <- "no event"
  n_t  <- n_f <- n_j <- 0


  #
  # Simulation

  j <- 1 # iteration
  repeat{

    # Compute next removal time
    not_recovered <- tau_J > t[j]
    tau_J_next    <- min(tau_J[not_recovered])

    # Generate candidate infection time
    tau_T_candidate <- if(S[j] > 0) {
      mu_j <- if(gener)  beta * S[j]^(1 - b) * I[j]  else beta * S[j] * I[j]
      t[j] + stats::rexp(1, mu_j)            # candidate infection time
    } else if(S[j] == 0) {            # susceptible pop depleted
      Inf
    }

    # Event: removal or infection
    if(tau_T_candidate < tau_J_next) { # infection (SIR) / exposition (SEIR)

      I_tau_t_true <- c(I_tau_t_true, I[j]) # sanity check for f_log()

      n_t <- n_t + 1
      X[j + 1] <- "t"
      t[j + 1] <- tau_T_candidate
      S[j + 1] <- S[j] - 1
      if(type == "SEIR") {
        E[j + 1] <- E[j] + 1
        I[j + 1] <- I[j]
        #tau_T[n_t] <- tau_T_candidate # check index
      } else if(type == "SIR") {
        I[j + 1] <- I[j] + 1
        tau_T[I0 + n_t] <- tau_T_candidate
      }

      # simulate removal time of newly infected
      iota <- simulate_iota(1, iota_dist, gamma, shape, rate)

      tau_J[I0 + n_t] <- tau_T[I0 + n_t] + iota

    } else if(tau_J_next <= tau_T_candidate) { # removal

      n_j <- n_j + 1
      X[j + 1] <- "j"
      t[j + 1] <- tau_J_next
      S[j + 1] <- S[j]
      if(type == "SEIR")  E[j + 1] <- E[j]
      I[j + 1] <- I[j] - 1

      if(I[j + 1] == 0)  break # infectious population depleted
      #TODO: if SEIR, also need to have E[j + 1]] == 0

    } # end-if

    j <- j + 1

  } # end-repeat


  # Events observed before t_end
  t_obs             <- t <= t_end
  tau_T_obs         <- tau_T <= t_end
  tau_J_obs         <- tau_J <= t_end
  tau_F_obs         <- tau_F <= t_end
  tau_T[!tau_T_obs] <- Inf
  tau_J[!tau_J_obs] <- Inf
  tau_F[!tau_F_obs] <- Inf

  x <- list(
    compatible = TRUE, tau_T = tau_T, tau_J = tau_J, tau_F = tau_F
  )

  # MLE
  n_t_obs    <- sum(0 < tau_T & is.finite(tau_T))
  n_j_obs    <- sum(is.finite(tau_J))

  dt              <- diff(c(t[t_obs], t_end))
  integral_I_obs  <- sum(I[t_obs] * dt)
  integral_SI_obs <- if(gener) {
    sum(I[t_obs] * S[t_obs]^(1 - b) * dt)
  } else{
    sum(I[t_obs] * S[t_obs]         * dt)
  }

  beta_MLE  <- n_t_obs / integral_SI_obs
  gamma_MLE <- n_j_obs / integral_I_obs
  #rate_MLE  <- n_j_obs / (sum(SS_true[["iota_recovered"]]^nu_true) + sum(SS_true[["iota_infectious"]]^nu_true))
  #shape_MLE <- optimize(
  # nu_posterior,
  # interval = c(0.01, 1e1),
  # theta = theta_true, SS = SS_true,
  # maximum = TRUE
  # )[["maximum"]]
  # TODO: run optim(), with true values as initial values to obtain MLE

  # Output
  out <- list(
    x = x, t = t, X = X, S = S, E = E, I = I, t_end = t_end, I0 = I0, S0 = S0,
    beta_MLE = beta_MLE, gamma_MLE = gamma_MLE,
    n_t_obs = n_t_obs, n_j_obs = n_j_obs,
    integral_SI_obs = integral_SI_obs, integral_I_obs = integral_I_obs
  )

  return(out)

}
