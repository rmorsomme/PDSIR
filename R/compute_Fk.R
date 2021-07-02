
#' Compute the number of expositions per observation interval
#'
#' @param x vector of infection times
#' @param ts observation schedule
#' @param S0 size of initial susceptible population
#'
#' @return number of infections per observation interval
#' @export
#'
compute_Fk <- function(x, ts, S0) { # only for SEIR

  tau_F <- x[["tau_F"]]
  tau_F <- tau_F[is.finite(tau_F) & tau_F > 0]

  K     <- length(ts) - 1
  F_k   <- rep(NA, K)

  for(k in 1 : K) {
    F_k[k] <- sum(dplyr::between(tau_F, ts[k], ts[k + 1]))
  }

  return(F_k)

}
