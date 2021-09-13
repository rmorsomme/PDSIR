
#' Compute the sufficient statistics from the latent space
#'
#' These sufficient statistics are used to conduct inference on the parameters
#'
#' @inheritParams rprop_x
#'
#' @param return_SI logical; whether to return the trajectories of S and I
#'
#' @return sufficient statistics
#' @export
#'
suff_stat <- function(x, Y, gener, b, return_SI = FALSE) {

  stopifnot(is.logical(gener), is.numeric(b), b > 0, is.logical(return_SI))

  # Verify compatibility of x
  if(!x[["compatible"]])  return(list(compatible = FALSE))

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
  iota_removed    <- tau_J[recovered] - tau_T[recovered]
  iota_infectious <- t_end            - tau_T[infectious]
  if(any(iota_removed < 1e-12))  return(list(compatible = FALSE))

  # Event time
  tau_T     <- tau_T[infected_during]
  tau_J     <- tau_J[recovered]
  tau       <- c(tau_T, tau_J)
  order_tau <- order(tau)

  # Compute S(tau), I(tau), I(tau_T)
  chi     <- c(rep(TRUE, n_T), rep(FALSE, n_J))[order_tau] # type of event (TRUE: infection, FALSE: recovery) # Sanity check: all.equal(X=="t", chi)

  delta_I <- ifelse(chi,  1, -1) # change in I
  delta_S <- ifelse(chi, -1,  0) # change in S

  I_tau   <- c(I0, I0 + cumsum(delta_I))
  S_tau   <- c(S0, S0 + cumsum(delta_S))
  I_tau_T <- I_tau[c(  chi, FALSE)] # include `FALSE` to exclude last element which corresponds to I(t_end)
  S_tau_T <- S_tau[c(  chi, FALSE)]

  if(any(I_tau_T == 0))  return(list(compatible = FALSE)) # depletion of infectious particles

  # Compute integrals
  tau         <- tau[order_tau]
  dtau        <- diff(c(0, tau, t_end))

  integral_SI <- if(gener)  sum(dtau * I_tau * S_tau^(1 - b))
                 else       sum(dtau * I_tau * S_tau        )
  integral_I  <- sum(dtau * I_tau)

  # Output
  if(return_SI)  return(list(S = S_tau, I = I_tau, t = c(0, tau))) # for plotting trajectories

  SS <- list( # for MCMC
    compatible = TRUE,
    n_T = n_T, n_J = n_J,
    iota_removed = iota_removed, iota_infectious = iota_infectious,
    integral_SI = integral_SI, integral_I = integral_I,
    I_tau_T = I_tau_T, S_tau_T = S_tau_T
  )
  return(SS)

}
