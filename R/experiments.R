#' Experiment 1: Proof of concept
#'
#' @param S0 size of initial susceptible population
#' @param I0 size initial infectious population
#' @param theta parameters of the SIR process
#' @param t_end end of observation period
#' @param K number of observation intervals
#' @param N number of iterations of the Markov chain
#' @param thin thinning parameter for the Markov chain
#' @param rho proportion of the latent data updated each iteration
#' @param param c("bg", "bR"); parameterize the models in terms of (beta, gamma) or (beta, R0)
#' @param approx c("poisson", "ldp"); whether to approximate the distribution of the infection times with a poisson process or a linear death process
#' @param iota_dist c("exponential", "weibull"); distribution of the infection period
#' @param gener logical; whether to use the generalized SIR of Severo (1972)
#' @param b parameter of the generalized SIR
#' @param path directory in which to save the figure
#' @param theta_0_factor factors by which the true value of the parameters is multiplied to initialize the Markov chain
#' @param plot_id name file for the figure(s)
#' @param save_fig logical; whether to save the figures generated
#'
#' @return list containing the parameters, observed data, Markov chain and run time of the algorithm
#' @export
#'
experiment_1_proof_of_concept <- function(
  S0 = 1e3, I0 = 1e1, theta = list(R0 = 2.5, lambda = 1, shape = 1, gamma = 1),
  t_end = 6, K = 10,
  N = 1e3, thin = 1, rho = 1/3,
  param = "bR", approx = "ldp",
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  theta_0_factor = 1,
  path = NULL, plot_id, save_fig = TRUE
  ) {

  theta <- complete_theta(theta, iota_dist, S0)

  # Observed data
  SIR   <- simulate_SEM(S0, I0, t_end, theta, iota_dist, gener, b)
  if(save_fig)  draw_trajectories(SIR, plot_id, path, t_end, type = "SIR")
  Y     <- observed_data(SIR, K)

  # run a long chain
  theta_0 <- purrr::map2(theta, theta_0_factor, `*`)
  if(iota_dist == "weibull")  theta_0[["shape"]] <- theta[["shape"]]

  MC <- run_DAMCMC(
    Y, N,
    rho, param, approx,
    iota_dist, gener, b,
    thin, theta_0
    )

  return(list(theta = theta, Y = Y, MC = MC, SIR = SIR))

}



#' Analyze output of MCMC run
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param x object returned by the function experiment_1_proof_of_concept
#' @param burnin number of iterations to discard, the default value NULL will discard half of the draws.
#' @param theta_true true value of the parameters
#' @param n_max last draw to take into account
#'
#' @return list of summary statistics with and without burn-in
#' @export
#'
experiment_1_output_analysis <- function(
  x, iota_dist = "exponential", theta_true,
  burnin = 0, thin = 1, n_max = NULL,
  plot_id = NULL, path = NULL, save_fig = TRUE
  ) {

  theta <- x[["theta"]]
  MC    <- x[["MC"   ]]
  Y     <- x[["Y"    ]]

  summary_no_burn <- analyze_MCMC(
    MC, burnin = 0, thin, n_max, iota_dist,
    plot_id = paste0(plot_id, "_no_burn"), path = path,
    theta_true = theta, Y = Y, save_fig = save_fig
    )

  if(is.null(burnin))  burnin <- length(MC[["theta"]])/2
  summary_burn <- analyze_MCMC(
    MC, burnin, thin, n_max, iota_dist,
    plot_id = paste0(plot_id, "_burn"), path = path,
    theta_true = theta, Y = Y, save_fig = save_fig
    )

  out <- list(summary_no_burn = summary_no_burn, summary_burn = summary_burn)
  return(out)

}


#' Experiment 2: Compare trajectories of SIR and PD-SIR process
#'
#' @inheritParams experiment_1_proof_of_concept
#' @inheritParams experiment_1_output_analysis
#'
#' @param Ks vector of numbers of time intervals to consider for the PD-SIR process
#'
#' @return figures comparing the trajectories of the PDSIR and SIR for different value K
#' @export
#'
experiment_2_PDSIR_trajectories <- function(
  S0 = 1e4, I0 = 1e1, theta = list(R0 = 4, gamma = 1),
  t_end = 6, Ks = c(3, 5, 10, 50, 1e3),
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  path
  ) {

  theta <- add_beta(theta, S0)

  # SIR
  SIR   <- simulate_SEM(S0, I0, t_end, theta)
  draw_trajectories(SIR, plot_id = "E2", path, t_end, type = "SIR")

  # PDSIR
  for(K in Ks) {

    Y        <- observed_data(SIR, K)
    PDSIR    <- rprop_x(theta, Y, gener, b, iota_dist, approx = "ldp")
    PDSIR_SI <- suff_stat(PDSIR, Y, gener, b, return_SI = TRUE)
    compare_trajectories(SIR, PDSIR_SI, paste0("E2_K", K), path, t_end)

  }

}



#' Experiment 3: Evaluate impact of rho on the M-H acceptance rate
#'
#' @inheritParams experiment_2_PDSIR_trajectories
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param S0s vector of sizes of initial susceptible population
#' @param R0s vector of R0's to be used in conjunction with S0s
#' @param gamma removal rate
#' @param rhos vector of values for the tuning parameter rho
#'
#' @return saves figures and a df summaryzing the results
#' @export
#'
experiment_3_acceptance_vs_rho <- function(
  S0s = c(1e2, 5e2, 1e3, 5e3), I0 = 1e1,
  R0s = c(2  , 2.5, 3  , 3.5), gamma = 1,
  t_end = 6, K = 20,
  rhos = c(0.02, 0.05, 0.1, 0.25, 0.5, 1),
  N = 1e3, thin = 1,
  iota_dist = "exponential",
  path
  ) {


  stopifnot(length(S0s) == length(R0s))

  results <- tibble::tibble(
    rho = numeric(), S0 = numeric(),
    accept_rate = numeric(), run_time = numeric(),
    ESS_beta = numeric(), ESS_gamma = numeric(), ESS_R0 = numeric(),
    ESSsec_beta = numeric(), ESSsec_gamma = numeric(), ESSsec_R0 = numeric()
    )

  for(k in 1 : length(S0s)) {

    # Parameters
    S0    <- S0s[k]
    theta <- list(R0 = R0s[k], gamma = gamma)

    # SEM
    SEM  <- simulate_SEM(S0, I0, t_end, theta)
    draw_trajectories(SEM, plot_id = paste0("E3_S0=", S0), path, t_end)
    Y    <- observed_data(SEM, K)

    for(rho in rhos) {

      MC <- run_DAMCMC(Y, N, rho, theta_0 = theta)
      print(paste0(S0, " - ", rho, ": ", Sys.time()))

      summary  <- analyze_MCMC(
        MC, burnin = min(N / 2, 1e4), iota_dist, thin,
        plot_id = paste0("E3_S0=", S0, "_rho=", rho), path,
        save_fig = FALSE,
        theta_true = theta,
        Y = Y
      )

      results  <- tibble::add_row(
        results,
        S0           = S0,
        rho          = rho,
        accept_rate  = summary[["rate_accept"]],
        run_time     = summary[["run_time"   ]],
        ESS_beta     = summary[["ESS"    ]][["beta" ]],
        ESS_gamma    = summary[["ESS"    ]][["gamma"]],
        ESS_R0       = summary[["ESS"    ]][["R0"   ]],
        ESSsec_beta  = summary[["ESS_sec"]][["beta" ]],
        ESSsec_gamma = summary[["ESS_sec"]][["gamma"]],
        ESSsec_R0    = summary[["ESS_sec"]][["R0"   ]]
      )

    } # end-for rho

  } # end-for S0

  # Output
  readr::write_csv(results, paste0(path, "/E3.csv"))

  # Figures
  draw_E3(results, path, "accept_rate" , "Acceptance Rate"  , "accept"     )
  draw_E3(results, path, "run_time"    , "Run Time"         , "runtime"    )
  draw_E3(results, path, "ESS_beta"    , "ESS for beta"     , "ESSbeta"    )
  draw_E3(results, path, "ESS_gamma"   , "ESS for gamma"    , "ESSgamma"   )
  draw_E3(results, path, "ESS_R0"      , "ESS for R0"       , "ESSR0"      )
  draw_E3(results, path, "ESSsec_beta" , "ESS/sec for beta" , "ESSsecbeta" )
  draw_E3(results, path, "ESSsec_gamma", "ESS/sec for gamma", "ESSsecgamma")
  draw_E3(results, path, "ESSsec_R0"   , "ESS/sec for R0"   , "ESSsecR0"   )

}

#' Experiment 4: estimate the coverage rate of the DA-MCMC algorithm
#'
#' @inheritParams experiment_3_acceptance_vs_rho
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param m number of independent iterations to estimate the coverage rate
#'
#' @return a tibble summarizing the coverage
#' @export
#'
experiment_4_coverage <- function(
  S0s = c(1e2, 5e2, 1e3), I0 = 1e1,
  R0s = c(2  , 2.5, 3  ), gamma = 1,
  t_end = 6, K = 10,
  N = 1e4, thin = 1, rho = 1,
  param = "bg", approx = "ldp",
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  theta_0_factor = 1,
  m = 1e2
) {

  stopifnot(length(S0s) == length(R0s))

  results <- tibble::tibble(
    S0 = numeric(), iteration = numeric(),
    var = numeric(), cover = numeric()
  )

  for(k in 1 : length(S0s)) {

    # Parameters
    S0    <- S0s[k]
    theta <- list(R0 = R0s[k], gamma = gamma)
    for(i in 1 : m) {
      print(paste0(S0, " - ", i, "/", m, ": ", Sys.time()))

      SEM <- simulate_SEM(S0, I0, t_end, theta)
      Y   <- observed_data(SEM, K)
      MC  <- run_DAMCMC(Y, N, rho, theta_0 = theta)

      summary  <- analyze_MCMC(
        MC, burnin = min(N / 2, 1e4), iota_dist, thin,
        plot_id = paste0("E3_S0=", S0, "_rho=", rho),
        save_fig = FALSE,
        theta_true = theta,
        Y = Y
      )

      results <- rbind(
        results,
        summary[["mean_quant_cover"]] %>%
          dplyr::select(.data$var, .data$cover) %>%
          dplyr::mutate(iteration = i, S0 = S0) %>%
          dplyr::select(S0, .data$iteration, .data$var, .data$cover)
      )

    } # end-for iteration
  } # end-for S0s

  return(results)

}


#' Experiment 5: Analysis of the Ebola data
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param theta_0 intial parameter values
#'
#' @return list containing the parameters, observed data, Markov chain and run time of the algorithm
#' @export
experiment_5_ebola <- function(
  I0 = 1e1, theta_0 = list(R0 = 1, gamma = 1e-1),
  N = 1e5, thin = 1, rho = 1,
  param = "bg", approx = "ldp",
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  path = NULL, plot_id = "E5"
) {

  #
  # Ebola Data

  d <- readr::read_csv(
    "../test/Data/Ebola/Guinea_modified.csv",
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
  T_k_ebola <- d[which(nam[[1]] == "GUECKEDOU"), ] # Ebola started in the GUECKEDOU prefecture

  # # Observed data (incidence)
  # print(T_k_ebola)
  # print(nam[[1]][14])
  # print(nam[[2]][c(1:3, 71:73)])
  #barplot(T_k_ebola, xlab = "Week")

  K     <- length(T_k_ebola) # number of time periods
  t_end <- 7 * K
  ts    <- seq(0, t_end, by = 7) # observation schedule

  # Observed data
  S0 <- 291823 # https://en.wikipedia.org/wiki/Prefectures_of_Guinea
  Y_ebola <- list(T_k = T_k_ebola, F_k = NULL, I0 = I0, S0 = S0, ts = ts, t_end = t_end)


  MC <- run_DAMCMC(
    Y_ebola, N,
    rho, param, approx,
    iota_dist, gener, b,
    thin, theta_0
  )

  summary <- analyze_MCMC(
    MC, burnin = min(1e4, N/2), iota_dist,
    plot_id = plot_id, path = path,
    theta_true = theta_0, Y = Y_ebola # TODO: add option to opt-out of theta_true
  )

  return(list(Y = Y_ebola, MC = MC, summary = summary))

}
