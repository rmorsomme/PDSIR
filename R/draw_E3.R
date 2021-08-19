
#' Line plots for experiment 3
#'
#' @inheritParams experiment_3_acceptance_vs_rho
#' @param results df of results from experiment 3
#' @param var variable of the df results to visualize
#' @param lab_y label of y-axis
#' @param plot_id name file for the figure
#'
#' @return saves a figure
#' @export

draw_E3 <- function(results, path, var, lab_y, plot_id) {

  results %>%
    ggplot2::ggplot(ggplot2::aes(.data$rho, .data[[var]], linetype = factor(.data$S0))) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(y = lab_y, x = expression(rho)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))

  ggplot2::ggsave(
    paste0("E3_", plot_id, ".jpg"),
    path = path, width = 1.61803, height = 1, scale = 5
    )

}
