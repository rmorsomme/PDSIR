draw_tp_hist_acf <- function(df, var, plot_id, path) {

  draw_traceplot(df, var, plot_id, path)
  draw_histogram(df, var, plot_id, path)
  draw_acf      (df, var, plot_id, path)

}
