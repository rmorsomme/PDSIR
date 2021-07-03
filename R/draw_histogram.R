
#' Generates a histogram
#'
#' @param df data frame
#' @param var variable under consideration
#' @param plot_id name file for the figure
#' @param bins number of bins in the histogram
#' @param path directory in which to save the figure
#'
#' @return a histogram
#' @export
#'
draw_histogram <- function(df, var, plot_id, path = NULL, bins = 30) {

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[var]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = .data$..count../sum(.data$..count..)), bins = bins) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))+
    ggplot2::labs(y = "Proportion")

  ggplot2::ggsave(
    paste(plot_id, var, "hist.jpg", sep = "_"),
    path = path, width = 1.61803, height = 1, scale = 5
  )

  return(g)

}
