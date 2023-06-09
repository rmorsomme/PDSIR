% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_SEM.R
\name{simulate_SEM}
\alias{simulate_SEM}
\title{Simulate trajectories of the Stochastic SIR and SEIR processes}
\usage{
simulate_SEM(
  S0 = 1000,
  I0 = 10,
  t_end = 6,
  theta = list(R0 = 2.5, gamma = 1, lambda = 1, shape = 1),
  iota_dist = "exponential",
  gener = FALSE,
  b = 1/2,
  E0 = 0,
  type = "SIR"
)
}
\arguments{
\item{S0}{size of initial susceptible population}

\item{I0}{size initial infectious population}

\item{t_end}{end of observation period}

\item{theta}{parameters of the SIR process}

\item{iota_dist}{c("exponential", "weibull"); distribution of the infection period}

\item{gener}{logical; whether to use the generalized SIR of Severo (1972)}

\item{b}{parameter of the generalized SIR}

\item{E0}{intial exposed population size}

\item{type}{c("SIR", "SEIR"); type of process to simulate}
}
\value{
a list with all types of useful objects
}
\description{
Also computes additional statistics such as MLEs and summary statistics
}
