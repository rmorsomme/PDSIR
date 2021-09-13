#' Draw trajectories of a SEM
#'
#' Line plot of the compartment sizes over time.
#'
#' @inheritParams analyze_MCMC
#' @inheritParams simulate_SEM
#'
#' @param SEM list corresponding to stochastic epidemic process
#'
#' @return save a figure of the compartments' trajectories
#' @export
#'

draw_trajectories <- function(
  SEM, plot_id = NULL, path, t_end = 10, type = "SIR"
) {

  t <- SEM[["t"]]
  S <- SEM[["S"]]
  E <- SEM[["E"]]
  I <- SEM[["I"]]
  X <- SEM[["X"]]
  N <- S[1] + I[1]

  tbl <- if(type == "SIR") {
    tibble::tibble(t = t, S = S,        I = I, R = N - S - I    , X = X)
  } else if(type == "SEIR") {
    tibble::tibble(t = t, S = S, E = E, I = I, R = N - S - E - I, X = X)
  }

  tbl %>%
    tidyr::pivot_longer(cols = (.data$S) : (.data$R), names_to = "Compartments", values_to = "count") %>%
    dplyr::filter(.data$t < t_end) %>%
    dplyr::mutate(Compartments = factor(.data$Compartments, levels = c("S", "E", "I", "R"))) %>%
    ggplot2::ggplot(ggplot2::aes(.data$t, .data$count, color = .data$Compartments)) +
    ggplot2::geom_line(size = 1.5) +
    ggplot2::theme(text = ggplot2::element_text(size = 25))

  ggplot2::ggsave(
    paste0(plot_id, "_trajectories.jpeg"),
    path = path, width = 1.5 * 1.61803, height = 1, scale = 5
  )

}
