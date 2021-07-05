
#' Simulate trajectories of the Stochastic SIR and SEIR processes
#'
#' Also computes additional statistics such as MLEs and summary statistics
#'
#' @param S0 size of initial susceptible population
#' @param I0 size initial infectious population
#' @param theta parameters of the SIR process
#' @param t_end end of observation period
#' @param generalized logical; whether to use generalized SIR
#' @param b parameter of generalized SIR
#' @param type c("SIR", "SEIR"); type of process to simulate
#'
#' @return list with all types of useful objects
#' @export
#'
simulate_SIR <- function(
  S0 = 1e3, I0 = 1e1, theta = list(R0 = 2, gamma = 1), t_end = 6,
  generalized = FALSE, b = 1/2,
  type = "SIR" # "SEIR"
) {


  # Setup
  N0 <- S0 + I0

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]
  beta  <- R0 / gamma / S0
  if(type == "SEIR")  epsilon <- theta[["epsilon"]]

  # Initialization
  tau_T <- tau_F <- tau_J <- rep(Inf, N0)
  S <- E <- I <- t <- W <- X <- I_tau_t_true <- c() # S, E, I, R, time, waiting time, type of event, and number of infectious at infection times

  if(type == "SEIR") {
    tau_J[1 : I0] <- NA
    tau_F[1 : I0] <- 0
  } else if(type == "SIR") {
    tau_J[1 : I0] <- 0
  }
  S[1] <- S0
  E[1] <- 0
  I[1] <- I0
  t[1] <- 0
  X[1] <- "no event"
  n_t  <- n_f <- n_j <- 0

  j <- 1 # iteration


  # Simulation
  repeat{

    # rates of events
    rate_T <- if(generalized)  beta * S[j]^(1 - b) * I[j]  else beta * S[j] * I[j]
    rate_F <- if(type == "SEIR")  epsilon * E[j]  else if(type == "SIR")  0
    rate_J <- gamma * I[j]
    rates  <- c(rate_T, rate_F, rate_J)

    # Draw time and type of next event
    W[j]     <- stats::rexp(1, sum(rates))       # waiting time until next event
    t[j + 1] <- t[j] + W[j]
    X[j + 1] <- sample(c("t", "f", "j"), 1, prob = rates) # type of event

    # update S, I and t (time)
    if(X[j + 1] == "t") { # infection / exposition (SEIR)

      n_t <- n_t + 1
      S[j + 1] <- S[j] - 1
      if(type == "SEIR") {
        E[j + 1] <- E[j] + 1
        I[j + 1] <- I[j]
        tau_T[n_t] <- t[j + 1]
      } else if(type == "SIR") {
        I[j + 1] <- I[j] + 1
        tau_T[I0 + n_t] <- t[j + 1]
      }
      I_tau_t_true <- c(I_tau_t_true, I[j]) # test: number of infectious at infection time


    } else if(X[j + 1] == "f") { # infection (SEIR)

      n_f <- n_f + 1
      S[j + 1] <- S[j]
      E[j + 1] <- E[j] - 1
      I[j + 1] <- I[j] + 1

      # Choose particle to be infected among exposed particles
      candidate <- which(is.finite(tau_T) & is.infinite(tau_F))
      infected  <- if(length(candidate) > 1)  sample(candidate, 1)  else  candidate
      tau_F[infected] <- t[j + 1]

    } else if(X[j + 1] == "j") { # removal

      n_j <- n_j + 1
      S[j + 1] <- S[j]
      I[j + 1] <- I[j] - 1
      if(type == "SEIR") {
        E[j + 1] <- E[j]
        candidate <- which(is.finite(tau_F) & is.infinite(tau_J))
      } else if(type == "SIR") {
        candidate <- which(is.finite(tau_T) & is.infinite(tau_J))
      }
      removed <- if(length(candidate) > 1)  sample(candidate, 1)  else  candidate
      tau_J[removed] <- t[j + 1]

      if(I[j + 1] == 0)  break  # no infectious remaining. Disease is over.

    } # end-if-else

    j <- j + 1
  } # end-repeat


  # Events observed before t_end
  t_obs             <- t <= t_end
  tau_T_obs         <- tau_T <= t_end
  tau_F_obs         <- tau_F <= t_end
  tau_J_obs         <- tau_J <= t_end
  tau_T[!tau_T_obs] <- Inf
  tau_F[!tau_F_obs] <- Inf
  tau_J[!tau_J_obs] <- Inf

  x <- list(
    compatible = TRUE, tau_T = tau_T, tau_J = tau_J, tau_F = tau_F
  )

  # MLE
  n_t_obs    <- sum(0 < tau_T & is.finite(tau_T))
  n_j_obs    <- sum(is.finite(tau_J))

  dt              <- diff(c(t[t_obs], t_end))
  integral_I_obs  <- sum(I[t_obs]            * dt)
  integral_SI_obs <- if(generalized)  sum(I[t_obs] * S[t_obs]^(1 - b) * dt)
  else                                sum(I[t_obs] * S[t_obs]         * dt)

  beta_MLE  <- n_t_obs / integral_SI_obs
  gamma_MLE <- n_j_obs / integral_I_obs

  # Output
  out <- list(
    x = x, t = t, X = X, S = S, E = E, I = I, t_end = t_end, I0 = I0, S0 = S0,
    beta_MLE = beta_MLE, gamma_MLE = gamma_MLE,
    n_t_obs = n_t_obs, n_j_obs = n_j_obs,
    integral_SI_obs = integral_SI_obs, integral_I_obs = integral_I_obs
  )

  return(out)

}
