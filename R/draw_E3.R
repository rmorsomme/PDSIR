
#' Line plot of M-H acceptance rate
#'
#' @inheritParams experiment_3_acceptance_vs_rho
#' @param results df of results from experiment 3
#'
#' @return saves a figure
#' @export
draw_accept_rate <- function(results, path) {

  results %>%
    ggplot2::ggplot(ggplot2::aes(.data$rho, .data$accept_rate, linetype = factor(.data$S0))) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(y = "Acceptance Rate",x = expression(rho)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))

  ggplot2::ggsave("E3_accept.jpg", path = path, width = 1.61803, height = 1, scale = 5)

}

draw_run_time <- function(results, path) {

  results %>%
    ggplot2::ggplot(ggplot2::aes(.data$rho, .data$run_time, linetype = factor(.data$S0))) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Run time (seconds)", x = expression(rho)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))

  ggplot2::ggsave("E3_runtime.jpg", path = path, width = 1.61803, height = 1, scale = 5)

}

draw_ESS_beta <- function(results, path) {

  results %>%
    ggplot(aes(rho, ESS_beta, linetype = factor(S0))) +
    geom_point() +
    geom_line() +
    labs(y = "Effective Sample Size", x = expression(rho)) +
    expand_limits(y = 0) +
    theme(text = element_text(size = 20))

  ggsave(paste0("E3_ESS_beta.jpg"), path = path, width = 1.61803, height = 1, scale = 5)

}


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
