
#' Compute the sufficient statistics from the latent space
#'
#' These sufficient statistics are used to conduct inference on the parameters
#'
#' @param x latent data
#' @param Y observed data
#' @param generalized logical; whether to use the generalized SIR
#' @param b parameter of the generalized SIR
#' @param return_SI logical; whether to return the trajectories of S and I
#'
#' @return sufficient statistics
#' @export
#'
sufficient_statistics <- function(x, Y, generalized, b, return_SI = FALSE) {

  # Verify compatibility of x
  if(!x[["compatible"]]) {
    return(list(compatible = FALSE))
  }

  # Setup
  tau_T <- x[["tau_T"]]
  tau_J <- x[["tau_J"]]
  t_end <- Y[["t_end"]]
  I0    <- Y[["I0"   ]]
  S0    <- Y[["S0"   ]]


  # Infected, infectious and recovered
  infected        <- is.finite(tau_T)
  infected_during <- infected & tau_T > 0 # excludes initially infectious
  recovered       <- is.finite(tau_J)
  infectious      <- infected & (! recovered)

  # Number of events
  n_T <- sum(infected_during)
  n_J <- sum(recovered)

  # Iota
  iota_recovered  <- tau_J[recovered] - tau_T[recovered]
  iota_infectious <- t_end            - tau_T[infectious]
  if(any(iota_recovered < 1e-12)) { # may happen that iota_recovered == 0 when rexp() produces 0
    return(list(compatible = FALSE))
  }

  # Event time
  tau_T     <- tau_T[infected_during]
  tau_J     <- tau_J[recovered]
  tau       <- c(tau_T, tau_J)
  order_tau <- order(tau)

  # S(tau), I(tau), I(tau_T)
  chi     <- c(rep(TRUE, n_T), rep(FALSE, n_J))[order_tau] # type of event (TRUE: infection, FALSE: recovery) # Sanity check: all.equal(X=="t", chi)

  delta_I <- ifelse(chi,  1, -1) # change in I
  delta_S <- ifelse(chi, -1,  0) # change in S

  I_tau   <- c(I0, I0 + cumsum(delta_I)) # Sanity check: all.equal(I, I_tau)
  S_tau   <- c(S0, S0 + cumsum(delta_S))
  I_tau_T <- I_tau[c(  chi, FALSE)] # include `FALSE` to exclude last element which corresponds to I(t_end) # all.equal(I_tau_t_true, I_tau_T)
  S_tau_T <- S_tau[c(  chi, FALSE)] # include `FALSE` to exclude last element which corresponds to I(t_end) # all.equal(I_tau_t_true, I_tau_T)

  if(any(I_tau_T == 0)) { # incompatible data
    return(list(compatible = FALSE))
  }

  # Compute integrals
  tau         <- tau[order_tau]         # times of event in order
  dtau        <- diff(c(0, tau, t_end)) # time between events # Sanity check: dtau - W

  integral_SI <- if(generalized)  sum(dtau * I_tau * S_tau^(1 - b)) # integral of step function
  else             sum(dtau * I_tau * S_tau        )
  integral_I  <- sum(dtau * I_tau)

  # Output
  if(return_SI) { # for plotting latent space
    return(list(S = S_tau, I = I_tau, t = c(0, tau)))
  }

  SS <- list( # for MCMC
    compatible = TRUE,
    n_T = n_T, n_J = n_J,
    iota_recovered = iota_recovered, iota_infectious = iota_infectious,
    integral_SI = integral_SI, integral_I = integral_I,
    I_tau_T = I_tau_T, S_tau_T = S_tau_T
  )
  return(SS)

}
