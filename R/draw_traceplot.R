
#' Draws a traceplot
#'
#' @inheritParams analyze_MCMC
#'
#' @param df data frame
#' @param var variable under consideration
#'
#' @return a traceplot
#' @export
#'
draw_traceplot <- function(df, var, plot_id, path = NULL) {

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$Iteration, y = .data[[var]])) +
    ggplot2::geom_line() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 40),
      axis.text.x = ggplot2::element_text(size = 30)
      )

  ggplot2::ggsave(
    paste(plot_id, var, "tp.jpg", sep = "_"),
    path = path, width = 1.61803, height = 1, scale = 5
  )

  return(g)

}
