
#

#' Removes the iterations of the transient phase of a Markov chain
#'
#' @param df data frame
#' @param b number of iterations to remove as burnin
#'
#' @return data frame without the burnin iterations
#' @export
#'
remove_burnin <- function(df, b = 0) {

  n <- nrow(df)
  df <- df %>%
    dplyr::filter(.data$Iteration > b, .data$Iteration < n) # second condition removes last iteration, which for LNA, is mistakenly set to 0

  return(df)

}
