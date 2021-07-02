

#' Probability of being removed before next observation time
#'
#' Compute p_i, the probability of being removed before t_end given an infection at time tau_T.
#'
#' @param theta named list with the parameters gamma (and nu and lambda)
#' @param tau_T observed event times
#' @param iota_distribution c("exponential", "weibull")
#' @param t_end time before which individuals may be removed
#'
#' @return a vector of the probabilities that individual infected at time tau_T are removed before t_end
#' @export
#'
compute_pi <- function(theta, tau_T, iota_distribution, t_end) {

  # Setup
  gamma  <- theta[["gamma"]]
  lambda <- theta[["lambda"]]
  nu     <- theta[["nu"]]

  # Compute p_i
  p_i <- if(iota_distribution == "exponential") {
    stats::pexp(t_end - tau_T, gamma)
  } else if(iota_distribution == "weibull") {
    pweibull2(t_end - tau_T, nu, rate = lambda)
  }

  # Output
  return(p_i)

}
