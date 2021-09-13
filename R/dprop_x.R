
#
#' Log proposal density for the latent data
#'
#' @inheritParams rprop_x
#'
#' @param i_update index set of particles whose values are updated
#'
#' @return log density
#' @export
#'
dprop_x <- function(
  theta, Y, x, i_update,
  gener, b, iota_dist, approx
) {

  # Setup

  T_k    <- Y[["T_k"  ]]
  I0     <- Y[["I0"   ]]
  ts     <- Y[["ts"   ]]
  t_end  <- Y[["t_end"]]
  K      <- length(T_k)
  T_sum  <- sum(T_k) + I0

  beta   <- theta[["beta"  ]]
  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  tau_T  <- x    [["tau_T"]]
  tau_J  <- x    [["tau_J"]]
  I_k    <- x    [["I_k"  ]]
  S_k    <- x    [["S_k"  ]]


  # Contribution of infections
  contribution_infection <- if(approx == "poisson") {
    0
  } else if(approx == "ldp") { # linear death process

    # exclude initially infectious particles
    i_update_tau_T <- setdiff(i_update, 1 : I0)

    mu_k  <- if(gener)  beta * I_k * S_k^(-b)  else  beta * I_k
    low   <- ts[1 : K      ] # TODO: compute only once and attach to Y
    upp   <- ts[2 : (K + 1)]

    mu_i  <- rep(mu_k, T_k)
    low_i <- rep(low , T_k)
    upp_i <- rep(upp , T_k)

    sum(dexp_trunc_log(
      tau_T[i_update_tau_T     ], mu_i [i_update_tau_T - I0],
      low_i[i_update_tau_T - I0], upp_i[i_update_tau_T - I0]
    ))
  } # end-if(approximation)

  # Contribution of removal
  p_i <- compute_pi(theta, tau_T, iota_dist, t_end)

  i_update_obs     <- setdiff(i_update, which(is.infinite(tau_J)))
  i_update_not_obs <- setdiff(i_update, which(is.finite  (tau_J)))

  contribution_removal <-
    sum(
      log(1 - p_i[i_update_not_obs])
    ) +
    sum(
      log(p_i[i_update_obs]) +
        contribution_observed_removal(
          theta, tau_T, tau_J, i_update_obs, iota_dist, t_end
        )
    )

  # Output
  loglik <- contribution_infection + contribution_removal
  return(loglik)

}
