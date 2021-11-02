#' Infection dynamics in the SIR and the PD-SIR processes.
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @return Figure contrasting the infection rate in the SIR and the PD-SIR processes.
#' @export
#'
draw_mu <- function(
  theta = list(R0 = 2, gamma = 10/2),
  S0 = 10, I0 = 2, t_end = 0.6, K = 3,
  iota_dist = "exponential", gener = FALSE, b = 1,
  path
  ) {

  set.seed(1)
  theta <- complete_theta(theta, iota_dist, S0)
  SIR   <- simulate_SEM(S0, I0, t_end, theta, iota_dist, gener, b)

  df_SIR <- tibble::tibble(
    t  = c(SIR$t, t_end),
    mu = c(theta$beta * SIR$I, 3) # add a row manually so that the steps do not stop with vertical line.
    )
  df_PD <- tibble::tibble(
    t     = seq(0, t_end, length.out = K + 1),
    t_SIR = purrr::map_dbl(.data$t, ~ max(df_SIR$t[df_SIR$t <= .]))
  ) %>%
    dplyr::left_join(df_SIR, by = c("t_SIR" = "t"))
  df_PD[K + 1, 3] <- df_PD[K, 3] # change last value of mu_PD so that step do not stop with vertical line.

  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$t, y = .data$mu)) +
    ggplot2::geom_step(data = df_SIR, alpha = 0.5, linetype = 2, size = 1) +
    ggplot2::geom_step(data = df_PD , alpha = 0.5, size = 1) +
    ggplot2::geom_vline(
      xintercept = SIR$t[SIR$X == "t"], col = "red", alpha = 0.4, linetype = 3, size = 1
    ) +
    ggplot2::geom_vline(
      xintercept = SIR$t[SIR$X == "j"], col = "blue", alpha = 0.4, linetype = 3, size = 1
    ) +
    ggplot2::xlim(0, t_end) +
    ggplot2::labs(y = "Infection rate")

  ggplot2::ggsave(
    "infection_rate_SIR_PDSIR.jpg",
    path = path, width = 1.61803, height = 1, scale = 4
  )

}
