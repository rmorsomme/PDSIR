% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propose_tau.R
\name{propose_tau_T}
\alias{propose_tau_T}
\title{Generates infection times}
\usage{
propose_tau_T(n, mu, lower, upper, approx)
}
\arguments{
\item{n}{number of infection times to generate}

\item{mu}{individual rate of infection}

\item{lower}{lower bound}

\item{upper}{upper bound}

\item{approx}{c("poisson", "ldp"); whether to approximate the distribution of the infection times with a poisson process or a linear death process}
}
\value{
infection times
}
\description{
Generates infection times
}
