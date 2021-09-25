
#' Line plots for experiment 3
#'
#' @inheritParams experiment_3_acceptance_vs_rho
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param results df of results from experiment 3
#' @param var variable of the df results to visualize
#' @param lab_y label of y-axis
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

#' Title
#'
#' @inheritParams draw_E3
#'
#' @return saves a figures
#' @export
#'
draw_E3_facet <- function(results, path, var, lab_y, plot_id) {

  results %>%
    ggplot2::ggplot(ggplot2::aes(.data$rho, .data[[var]])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(y = lab_y, x = expression(rho)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::facet_grid(S0 ~ ., labeller = ggplot2::label_both, scales = "free")

  ggplot2::ggsave(
    paste0("E3_facet_", plot_id, ".jpg"),
    path = path, width = 1.61803, height = 1, scale = 5
  )

}
