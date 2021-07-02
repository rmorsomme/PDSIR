
#' Discrete Incidence Data for Infections
#'
#' Compute the number of infections in each time interval from the output of a MCMC
#'
#' @param x vector of infection times
#' @param ts observation schedule
#' @param S0 initial population of susceptibles
#'
#' @return a vector of the number of infections in each time interval
#' @export
#'
#'
compute_Tk <- function(x, ts, S0) {

  tau_T <- x[["tau_T"]]
  tau_T <- tau_T[is.finite(tau_T) & tau_T > 0] # exclude infinite values and zeros

  K     <- length(ts) - 1
  T_k   <- rep(NA, K)

  for(k in 1 : K) {
    T_k[k] <- sum(dplyr::between(tau_T, ts[k], ts[k + 1]))
  }

  return(T_k)

}



#' Discrete Incidence Data for Infections (Fast)
#'
#' @param x vector of infection times
#' @param ts observation schedule
#' @param S0 initial population of susceptibles
#'
#' @return a vector of the number of infections in each time interval
#' @export
#'
compute_Tk_MH_fast <- function(x, ts, S0) {

  tau_T <- x[["tau_T"]]

  K     <- length(ts) - 1
  T_k   <- rep(NA, K)

  for(k in 1 : K) {
    T_k[k] <- sum(dplyr::between(tau_T, ts[k], ts[k + 1]))
  }

  return(T_k)

}
