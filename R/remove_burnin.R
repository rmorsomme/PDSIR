
#

#' Removes the iterations of the transient phase of a Markov chain
#'
#' @param df data frame
#' @param burnin number of iterations to remove as burnin
#'
#' @return data frame without the burnin iterations
#' @export
#'
remove_burnin <- function(df, burnin = 0) {

  df <- df %>%
    dplyr::filter(.data$Iteration > burnin)

  return(df)

}
