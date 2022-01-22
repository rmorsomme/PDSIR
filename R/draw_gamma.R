#' Illustrates minorization of gamma densities in manuscript
#'
#' @param path directory in which to save the figure
#'
#' @return Illustration of gamma minorization
#' @export
draw_gamma <- function(path){

  # Setup
  x <- seq(0, 8, by = 0.01)

  a <- 2
  b <- 0.5
  A <- 1
  B <- 1

  xa <- a/B*log(1+B/b)
  xb <- 1/b*(gamma(a+A)/gamma(a))^(1/A)
  xA <- (a+A)/B*log(1+B/b)
  xB <- 1/(b+B)*(gamma(a+A)/gamma(a))^(1/A)


  # Gamma 1

  df1 <- expand.grid(alpha = seq(0, A, by = 0.25), beta = 0, x = x) %>%
    dplyr::mutate(
      Density = list(.data$x, .data$alpha, .data$beta) %>%
        purrr::pmap(~ dgamma(..1, a + ..2, b + ..3))
      ) %>%
    tidyr::unnest(.data$Density)

  ggplot2::ggplot(df1, ggplot2::aes(x = .data$x, y = .data$Density)) +
    ggplot2::geom_line(ggplot2::aes(col = as.factor(.data$alpha)), size = 1) +
    ggplot2::scale_color_discrete(name = expression(Parameter~alpha)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::geom_vline(xintercept = xb) +
    ggplot2::annotate(geom = "text", x = xb + 0.35, y = 0.18, label = expression(x[b]), size = 6)

  ggplot2::ggsave(
    "gamma1.jpg",
    path = path, width = 1.61803, height = 1, scale = 4
  )


  # Gamma 2
  df2 <- expand.grid(alpha = 0, beta = seq(0, B, by = 0.25), x = x) %>%
    dplyr::mutate(
      Density = list(.data$x, .data$alpha, .data$beta)
      %>% purrr::pmap(~ dgamma(..1, a + ..2, b + ..3))) %>%
    tidyr::unnest(.data$Density)

  ggplot2::ggplot(df2, ggplot2::aes(x = .data$x, y = .data$Density)) +
    ggplot2::geom_line(ggplot2::aes(col = as.factor(beta)), size = 1) +
    ggplot2::scale_color_discrete(name = expression(Parameter~beta)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::geom_vline(xintercept = xa) +
    ggplot2::annotate(geom = "text", x = xa + 0.35, y = 0.525, label = expression(x[a]), size = 6)

  ggplot2::ggsave(
    "gamma2.jpg",
    path = path, width = 1.61803, height = 1, scale = 4
  )

  # Gamma joint
  df3 <- expand.grid(alpha = c(0,1), beta = c(0,1), x = x) %>%
    dplyr::mutate(
      Density = list(.data$x, .data$alpha, .data$beta) %>%
        purrr::pmap(~ dgamma(..1, a + ..2, b + ..3))) %>%
    tidyr::unnest(.data$Density) %>%
    dplyr::mutate(code = paste0(.data$alpha, .data$beta))

  ggplot2::ggplot(df3, ggplot2::aes(x = .data$x, y = .data$Density)) +
    ggplot2::geom_line(ggplot2::aes(col = .data$code), size = 1) +
    ggplot2::scale_color_discrete(
      name = expression(paste("Parameters (", .data$alpha, ", " ,.data$beta, ")")),
      breaks = c("00", "01", "10", "11"),
      labels = c("(0, 0)", "(0, 1)", "(1, 0)", "(1, 1)")
      ) +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::geom_vline(xintercept = c(xa, xA, xb, xB)) +
    ggplot2::annotate(geom = "text", x = xa + 0.35, y = 0.525, label = expression(x[a]), size = 6) +
    ggplot2::annotate(geom = "text", x = xb + 0.35, y = 0.525, label = expression(x[b]), size = 6) +
    ggplot2::annotate(geom = "text", x = xA + 0.35, y = 0.525, label = expression(x[A]), size = 6) +
    ggplot2::annotate(geom = "text", x = xB + 0.35, y = 0.525, label = expression(x[B]), size = 6)

  ggplot2::ggsave(
    "gamma3.jpg",
    path = path, width = 1.61803, height = 1, scale = 4
  )

}
