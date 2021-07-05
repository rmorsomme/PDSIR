
#' Generates a bi-dimensional density plot
#'
#' @inheritParams draw_traceplot
#'
#' @param var1 variable to plot on the x-axis
#' @param var2 variable to plot on the y-axis
#' @param xlim vector of length 2; limits on the x-axis
#'
#' @return a figure of a biavriate density
#' @export
#'
draw_density_2d <- function(df, var1, var2, plot_id, path = NULL, xlim = c(NA_real_, NA_real_)) {

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[var1]], y = .data[[var2]])) +
    ggplot2::geom_density_2d_filled() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::xlim(xlim[[1]], xlim[[2]])

  ggplot2::ggsave(
    paste(
      plot_id, var1, var2, "bi_dens.jpg",
      sep = "_"
      ),
    path = path, width = 1.61803, height = 1, scale = 5
  )

  return(g)
}
