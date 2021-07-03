
#' Creates figures and summary statistics of the output of the Markov chain
#'
#' @param MC output of the MCMC algorithm
#' @param burnin number of iterations to discard
#' @param thin thinning argument for the iterations of the Markov chain
#' @param plot_id name file for the figures
#' @param save_fig logical; whether to save the figures generated
#' @param do_SS logical; whether to analyze the summary statistics in addition to the parameters
#' @param beta_true true value of the beta parameter
#' @param gamma_true true value of the gamma parameter
#' @param path directory in which to save the figure
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return a list of summary statistics including posterior means, posterior quantiles, effective sample size and the acceptance rate
#' @export
#'
analyze_MCMC <- function(
  MC, burnin = 0, thin = 1,
  plot_id = NULL, path = NULL, save_fig = TRUE, do_SS = FALSE,
  beta_true, gamma_true
) {

  # Setup
  theta         <- MC[["theta"      ]]
  loglik        <- MC[["loglik"     ]]
  if(do_SS)  SS <- MC[["SS"         ]]
  rate_accept   <- MC[["rate_accept"]]
  S0            <- MC[["S0"         ]]

  # Data Wrangling
  theta_tidy <- data.table::rbindlist(theta) %>%
    add_iteration %>%
    remove_burnin(burnin) %>%
    dplyr::filter(.data$Iteration %% thin == 0)

  if(do_SS) {
    SS_tidy     <- SS %>%
      purrr::map( ~ .[c("n_T", "n_J", "integral_SI", "integral_I")]) %>%
      data.table::rbindlist  %>%
      add_iteration %>%
      remove_burnin(burnin)  %>%
      dplyr::filter(.data$Iteration %% thin == 0)
  }

  # Figures
  if(save_fig) {
    # Traceplots
    tp_beta  <- draw_traceplot(theta_tidy, "beta"       , plot_id, path)
    tp_gamma <- draw_traceplot(theta_tidy, "gamma"      , plot_id, path)
    tp_R0    <- draw_traceplot(theta_tidy, "R0"         , plot_id, path)
    if(do_SS) {
      #tp_nt    <- draw_traceplot(SS_tidy   , "n_T"        , plot_id, path)
      tp_nj    <- draw_traceplot(SS_tidy   , "n_J"        , plot_id, path)
      tp_I     <- draw_traceplot(SS_tidy   , "integral_I" , plot_id, path)
      tp_SI    <- draw_traceplot(SS_tidy   , "integral_SI", plot_id, path)
    }
    tp_ll    <- tibble::tibble(loglik = loglik) %>%
      add_iteration %>%
      remove_burnin(burnin) %>%
      draw_traceplot("loglik", plot_id, path)

    # Histograms
    hist_beta  <- draw_histogram(theta_tidy, "beta"       , plot_id, path)
    hist_gamma <- draw_histogram(theta_tidy, "gamma"      , plot_id, path)
    hist_R0    <- draw_histogram(theta_tidy, "R0"         , plot_id, path)
    if(do_SS) {
      hist_nt    <- draw_histogram(SS_tidy   , "n_T"        , plot_id, path)
      hist_nj    <- draw_histogram(SS_tidy   , "n_J"        , plot_id, path)
      hist_I     <- draw_histogram(SS_tidy   , "integral_I" , plot_id, path)
      hist_SI    <- draw_histogram(SS_tidy   , "integral_SI", plot_id, path)
    }

    # Joint Densities
    den_bg     <- draw_density_2d(theta_tidy, "beta", "gamma", plot_id, path)
    den_Rg     <- draw_density_2d(theta_tidy, "R0"  , "gamma", plot_id, path)
    den_Rb     <- draw_density_2d(theta_tidy, "R0"  , "beta" , plot_id, path)
    # draw_hpd_2d(theta_tidy, "beta", "gamma", beta_true, gamma_true, plot_id) # too buggy, need to find a better alternative

    # ACF
    acf_beta  <- draw_acf(theta_tidy, "beta" , plot_id, path = path)
    acf_gamma <- draw_acf(theta_tidy, "gamma", plot_id, path = path)
    acf_R0    <- draw_acf(theta_tidy, "R0"   , plot_id, path = path)

  }

  # Summary Statistics
  out <- list(
    post_mean   = colMeans(dplyr::select(theta_tidy, - .data$Iteration)),
    quant_05    = dplyr::select(theta_tidy, - .data$Iteration) %>% purrr::map_dbl(stats::quantile, probs = 0.05),
    quant_95    = dplyr::select(theta_tidy, - .data$Iteration) %>% purrr::map_dbl(stats::quantile, probs = 0.95),
    rate_accept = rate_accept,
    ESS         = coda::effectiveSize(coda::mcmc(theta_tidy))
  )

  return(out)

}
