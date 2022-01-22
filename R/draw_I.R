
#' Barplot of the Observed Infection Counts
#'
#' @param path_data directory in which the .csv file is located
#' @param path_figure directory in which to save the figure
#'
#' @return Barplot of infection counts.
#' @export
#'
draw_Ebola_data <- function(
  path_data = "Input/Ebola/Guinea_modified.csv",
  path_figure
  ) {

  data <- extract_Ebola_data(path_data)
  I_K <- data[["I_K"]]
  ts  <- data[["ts"]]
  ts <- ts[-1]

  tibble::tibble(I = I_K, t = ts) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$t, y = .data$I)) +
    ggplot2::geom_bar()

  ggplot2::ggsave(
    "infection_rate_SIR_PDSIR.jpg",
    path = path_figure, width = 1.61803, height = 1, scale = 5
  )

}

#'
#' Extract Observed Infection Counts from .csv File
#'
#' @param path directory in which the .csv file is located
#'
#' @return Infection counts for Gueckedou
#' @export
#'
extract_Ebola_data <- function(path) {

  d <- readr::read_csv(
    path,
    col_types = readr::cols(
      .default            = readr::col_double(),
      Location            = readr::col_character(),
      `Ebola data source` = readr::col_character(),
      `Indicator type`    = readr::col_character(),
      `Case definition`   = readr::col_character()
    )
  ) %>%
    dplyr::select(-(2:4)) %>%
    dplyr::group_by(.data$Location) %>%
    dplyr::summarize_all(sum) %>%
    tibble::column_to_rownames("Location") %>%
    as.matrix

  nam <- dimnames(d)
  dimnames(d) <- NULL
  #d[, 1 : 10]
  I_k_ebola <- d[which(nam[[1]] == "GUECKEDOU"), ] # Ebola started in the GUECKEDOU prefecture


  K     <- length(I_k_ebola) # number of time periods
  t_end <- 7 * K
  ts    <- seq(0, t_end, by = 7) # observation schedule

  return(list(I_K = I_k_ebola, ts = ts))

}
