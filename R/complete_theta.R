#' Include beta to the parameter vector theta
#'
#' @param theta parameters (without beta)
#' @param S0 initial population size
#'
#' @return theta including beta
#' @export
#'
add_beta <- function(theta, S0) {

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]

  theta[["beta"]] <- R0 * gamma / S0

  return(theta)

}

add_gamma <- function(theta) {

  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  theta[["gamma"]] <- 1 / (lambda ^ (- 1 / shape) * gamma(1 + 1 / shape))

  return(theta)

}


complete_theta <- function(theta, iota_dist, S0) {

  if(iota_dist == "exponential") theta[c("shape", "lambda")] <- NULL
  if(iota_dist == "weibull")     theta <- add_gamma(theta)
  theta <- add_beta(theta, S0)

  return(theta)

}
