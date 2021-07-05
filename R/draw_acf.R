
#' Generate a plot of the auto-correlation function
#'
#' @inheritParams draw_traceplot
#'
#' @param lag_max largest lag shown
#'
#' @return Barplot of the ACF
#' @export
#'
draw_acf <- function(df, var, plot_id, path = NULL, lag_max = 1.5e3) {

  x     <- dplyr::pull(df, .data[[var]])
  x_acf <- stats::acf(x, lag.max = lag_max, plot = FALSE)
  df    <- tibble::tibble(ACF = as.numeric(x_acf[["acf"]]), Lag = as.numeric(x_acf[["lag"]]))

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(.data$Lag, .data$ACF)) +
    ggplot2::geom_col(width = .25) +
    ggplot2::theme(text = ggplot2::element_text(size = 30))

  ggplot2::ggsave(
    paste(plot_id, var, "acf.jpg", sep = "_"),
    path = path, width = 1.61803, height = 1, scale = 5
  )

  return(g)

}
