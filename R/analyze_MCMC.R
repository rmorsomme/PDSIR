
#' Creates figures and summary statistics of the output of the Markov chain
#'
#' @inheritParams run_DAMCMC
#' @inheritParams experiment_1_output_analysis
#'
#' @param MC output of the MCMC algorithm
#' @param do_SS logical; whether to analyze the summary statistics in addition to the parameters
#' @param theta_true true value of the parameters
#' @param coverage whether to determine if the credible intervals contain the true value of the parameters
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @return a list of summary statistics including posterior means, posterior quantiles, effective sample size and the acceptance rate
#' @export
#'
analyze_MCMC <- function(
  MC, burnin = 0, thin, n_max = NULL, iota_dist,
  plot_id = NULL, path = NULL, save_fig = TRUE, do_SS = FALSE,
  theta_true, Y, coverage = TRUE
  ) {


  # Setup
  theta         <- MC[["theta"      ]]
  loglik        <- MC[["loglik"     ]]
  if(do_SS)  SS <- MC[["SS"         ]]
  rate_accept   <- MC[["rate_accept"]]
  S0            <- MC[["S0"         ]]
  run_time      <- MC[["run_time"   ]]

  S0            <- Y [["S0"         ]]

  if(is.null(n_max)) n_max <- length(theta) * thin

  # Data Wrangling
  theta_tidy <- data.table::rbindlist(theta) %>%
    add_iteration(thin) %>%
    dplyr::mutate(loglik = loglik) %>%
    dplyr::filter(.data$Iteration <= n_max) %>%
    remove_burnin(burnin) %>%
    dplyr::mutate(expected_infection_length = 1 / gamma)

  if(do_SS) {
    SS_tidy <- SS %>%
      purrr::map( ~ .[c("n_T", "n_J", "integral_SI", "integral_I")]) %>%
      data.table::rbindlist %>%
      add_iteration %>%
      remove_burnin(burnin)  %>%
      dplyr::filter(.data$Iteration %% thin == 0)
  } # end-if do_SS

  # Figures
  if(save_fig) {

    # Traceplots, histograms, ACF
    vars <- c("beta", "gamma", "R0", "loglik", "expected_infection_length")
    if(iota_dist == "weibull") vars <- c(vars, "lambda")
    if(do_SS) vars <- c(vars, c("n_J", "integral_I", "integral_SI"))

    for(var in vars)
      draw_tp_hist_acf(theta_tidy, var, plot_id, path)

    # Joint Densities
    vars_x <- c("R0"  , "R0"   , "beta" )
    vars_y <- c("beta", "gamma", "gamma")
    if(iota_dist == "weibull"){
      vars_x <- c(vars_x, "beta"  , "R0"    )
      vars_y <- c(vars_y, "lambda", "lambda")
    }

    for(i in 1 : length(vars_x))
      draw_density_2d(theta_tidy, vars_x[i], vars_y[i], plot_id, path)

    # draw_hpd_2d(theta_tidy, "beta", "gamma", beta_true, gamma_true, plot_id) # too buggy, need to find a better alternative

  } # end-if save_fig

  # Posterior mean
  post_mean <- theta_tidy %>%
    dplyr::select(- .data$Iteration, - .data$loglik) %>%
    {tibble::tibble(
      var  = colnames(.),
      mean = colMeans(.)
    )}

  # Output
  out <- list(
    run_time    = run_time,
    post_stat   = post_mean,
    rate_accept = rate_accept,
    ESS         = coda::effectiveSize(coda::mcmc(theta_tidy))
  )

  out[["ESS_sec"]] <- out[["ESS"]] / run_time


  # Coverage
  if(coverage){

    theta_true_df   <- theta_true %>%
      complete_theta(iota_dist, S0) %>%
      {tibble::tibble(var = names(.), theta_true = as.double(.))}

    post_mean_quant_cover <- theta_tidy %>%
      dplyr::select(- .data$Iteration, - .data$loglik) %>%
      {tibble::tibble(
        var  = colnames(.),
        mean = colMeans(.),
        quant_low = purrr::map_dbl(., stats::quantile, probs = 0.05),
        quant_upp = purrr::map_dbl(., stats::quantile, probs = 0.95)
      )} %>%
      dplyr::left_join(theta_true_df, by = "var") %>%
      dplyr::mutate(
        cover = list(.data$theta_true, .data$quant_low, .data$quant_upp) %>%
          purrr::pmap_lgl(dplyr::between)
      )

    out[["post_stat"]] <- post_mean_quant_cover

  }

  return(out)

}
