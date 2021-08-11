
#' Re-parameterization for the weibull distribution
#'
#' Provides the scale parameter of the Weibull distribution
#'
#' @inheritParams pweibull2
#'
#' @return scale parameter
#' @export
#'
rate_to_scale <- function(rate, shape){
  scale <- rate ^ (- 1 / shape)
  return(scale)
}
