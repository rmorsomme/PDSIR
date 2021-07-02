#' Title
#'
#' @param x vector of one or more numeric values
#' @param low lower bound
#' @param upp upper bound
#'
#' @return logical vector the same length as x
#' @export
#'
between_hard <- function(x, low, upp) {

  stopifnot(low <= upp, "low must be smaller than upp")

  is_between <- (low < x) & (x < upp)
  return(is_between)

}
