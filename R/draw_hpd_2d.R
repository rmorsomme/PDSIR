
#' draws highest posterior density
#' need to find better alternative because the HPDregionplot function is very buggy. check Peter Hoff's (2009) book
#'
#' @param df data frame
#' @param var1 variable to plot on the x-axis
#' @param var2 variable to plot on the y-axis
#' @param true1 true value of the parameter corresponding to variable 1
#' @param true2 true value of the parameter corresponding to variable 2
#' @param plot_id name file for the figures
#' @param path directory in which to save the figure
#' @param levels c(0.50, 0.95); level(s) at which to draw the HPD regions
#'
#' @return Contour plot with HPD region(s) at the specified level(s)
#' @export
#'
draw_hpd_2d <- function(df, var1, var2, true1, true2, plot_id, path = NULL, levels) {

  df_mcmc <- df %>%
    dplyr::select(.data[[var1]], .data[[var2]]) %>%
    coda::as.mcmc()

  # Figure
  grDevices::jpeg(file = paste0(path, "/", plot_id, "_HPD.jpeg"))
  emdbook::HPDregionplot(
    df_mcmc, n = 500,
    prob       = levels,
    xlab       = expression(beta), ylab = expression(gamma),
    cex.lab    = 2
  )
  graphics::points(true1, true2, pch = 8, col = "red")
  grDevices::dev.off()

}
