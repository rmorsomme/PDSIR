#' PDSIR: A package for fitting the stochastic SIR process with incidence data
#'   with a data-augmentation MCMC algorithm that uses the piece-wise decoupled SIR process
#'   to proposes the latent data.
#'
#' The PDSIR package provides all the functions necessary to run the DA-MCMC algorithm
#'   and analyze the output.
#'
#' @section PD-SIR:
#' The function \code{rprop_x} generates a process from the PD-SIR
#'  which by construction is compatible with the observed incidence data.
#'
#' @section DA-MCMC:
#' The function \code{run_DAMCMC} runs the DA-MCMC algorithm.
#'
#' @docType package
#' @name PDSIR
NULL
