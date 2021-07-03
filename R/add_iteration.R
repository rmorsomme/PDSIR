

#' Add an "Iteration" column to a data frame
#'
#' @param df data frame
#'
#' @return data frame with a variable for iteration
#' @export
#'
add_iteration <- function(df) {

  n  <- nrow(df)
  df <- dplyr::mutate(df, Iteration = 1 : n)
  return(df)

}
