

#' Add an "Iteration" column to a data frame
#'
#' @inheritParams analyze_MCMC
#'
#' @param df data frame
#'
#' @return data frame with a variable for iteration
#' @export
#'
add_iteration <- function(df, thin) {


  stopifnot(is.data.frame(df))

  n  <- nrow(df)
  df <- dplyr::mutate(df, Iteration = thin * (1:n))
  return(df)

}
