
#' Draw trajectories of the PDSIR and SIR
#'
#' Compare the compartment trajectories of the stochastic SIR process and of the PD-SIR process
#'
#' @inheritParams experiment_2_PDSIR_trajectories
#'
#' @param SIR_SI compartment sizes of SIR process
#' @param PDSIR_SI compartment sizes of PD-SIR process
#' @param plot_id name file for the figures
#'
#' @return save figures of trajectories
#' @export
#'
compare_trajectories <- function(SIR_SI, PDSIR_SI, plot_id = NULL, path, t_end){

  N <- SIR_SI[["S"]][1] + SIR_SI[["I"]][1]

  SIR <- tibble::tibble(
    t = SIR_SI[["t"]],
    S = SIR_SI[["S"]],
    I = SIR_SI[["I"]],
    R = N - .data$S - .data$I
      ) %>%
    tidyr::pivot_longer(
      cols = .data$S : .data$R, names_to = "Compartments", values_to = "count"
      ) %>%
    dplyr::filter(.data$t < t_end) %>%
    dplyr::mutate(Compartments = factor(.data$Compartments, levels = c("S", "I", "R")))

  PD_SIR <- tibble::tibble(
    t = PDSIR_SI[["t"]],
    S = PDSIR_SI[["S"]],
    I = PDSIR_SI[["I"]],
    R = N - .data$S - .data$I
      ) %>%
    tidyr::pivot_longer(
      cols = .data$S : .data$R, names_to = "Compartments", values_to = "count"
      ) %>%
    dplyr::filter(.data$t < t_end) %>%
    dplyr::mutate(Compartments = factor(.data$Compartments, levels = c("S", "I", "R")))

  ggplot2::ggplot(mapping = ggplot2::aes(
    .data$t, .data$count, group = .data$Compartments, color = .data$Compartments
    )) +
    ggplot2::geom_line(data = SIR, alpha = 0.25, color = "black", size = 1.25) +
    ggplot2::geom_line(data = PD_SIR, size = 1.25) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 30),
      legend.position = "none"
      )

  ggplot2::ggsave(
    paste0(plot_id, ".jpeg"),
    path = path, width = 1.5 * 1.61803, height = 1, scale = 5
  )

}
