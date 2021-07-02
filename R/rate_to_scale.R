
#

#' Re-parameterization for the weibull distribution
#'
#' Provides the scale parameter of the Weibull distribution
#'
#' @param rate rate parameter
#' @param shape shape parameter
#'
#' @return scale parameter
#' @export
#'
rate_to_scale <- function(rate, shape){
  scale <- rate ^ (- 1 / shape)
  return(scale)
}
