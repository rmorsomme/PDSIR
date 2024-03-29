% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEM-simulate.R
\name{observed_data}
\alias{observed_data}
\title{Compute observed data from trajectories of the stochastic SIR process}
\usage{
observed_data(SEM, K = 10, ts = NULL)
}
\arguments{
\item{SEM}{complete data for a stochastic epidemic model}

\item{K}{number of observation intervals}

\item{ts}{observation schedule}
}
\value{
list of observed data
}
\description{
Compute observed data from trajectories of the stochastic SIR process
}
