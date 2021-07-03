
#' Generates infection times
#'
#' @param n number of infection times to generate
#' @param mu individual rate of infection
#' @param lower lower bound
#' @param upper upper bound
#' @param approximation c("poisson", "pld"); whether to approximate the distribution of the infection times with a poisson process or a linear death process
#'
#' @return infection times
#' @export
#'
propose_tau_T <- function(n, mu, lower, upper, approximation) {

  tau_T <- if(approximation == "poisson") {
    stats::runif(n,     lower, upper)
  } else if(approximation == "pld") {
    rexp_trunc  (n, mu, lower, upper)
  }

  return(tau_T)

}
