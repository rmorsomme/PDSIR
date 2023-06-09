% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dprop_x.R
\name{dprop_x}
\alias{dprop_x}
\title{Log proposal density for the latent data}
\usage{
dprop_x(theta, Y, x, i_update, gener, b, iota_dist, approx)
}
\arguments{
\item{theta}{current values of the parameters for the SIR}

\item{Y}{observed data}

\item{x}{current configuration of the latent data}

\item{i_update}{index set of particles whose values are updated}

\item{gener}{logical; whether to use the generalized SIR of Severo (1972)}

\item{b}{parameter of the generalized SIR}

\item{iota_dist}{c("exponential", "weibull"); distribution of the infection period}

\item{approx}{c("poisson", "ldp"); whether to approximate the distribution of the infection times with a poisson process or a linear death process}
}
\value{
log density
}
\description{
Log proposal density for the latent data
}
