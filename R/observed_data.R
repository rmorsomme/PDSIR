
#' Compute observed data from trajectories of the stochastic SIR process
#'
#' @param SEM complete data for a stochastic epidemic model
#' @param K number of observation intervals
#' @param ts observation schedule
#' @param type c("SIR", "SEIR"); type of process to simulate
#'
#' @return list of observed data
#' @export
#'
observed_data <- function(SEM, K = 10, ts = NULL, type = "SIR") {

  # Setup
  t_end <- SEM[["t_end"]]
  x     <- SEM[["x"]]
  I0    <- SEM[["I0"]]
  S0    <- SEM[["S0"]]

  if(is.null(ts)){
    dt <- t_end / K
    ts <- seq(0, t_end, by = dt) # observation schedule (need not be regular)
  }

  # Observed data
  T_k <- compute_Tk(x, ts = ts             )
  if(type == "SEIR") F_k <- compute_Fk(x, ts = ts)

  # Output
  out <- list(T_k = T_k, F_k = F_k, I0 = I0, S0 = S0, ts = ts, t_end = t_end)
  return(out)

}
