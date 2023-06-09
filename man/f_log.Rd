% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_log.R
\name{f_log}
\alias{f_log}
\title{Log likelihood of the stochastic SIR model}
\usage{
f_log(theta, SS, gener, b, iota_dist = "exponential")
}
\arguments{
\item{theta}{current values of the parameters for the SIR}

\item{SS}{sufficient statistics of the current configuration of the latent data}

\item{gener}{logical; whether to use the generalized SIR of Severo (1972)}

\item{b}{parameter of the generalized SIR}

\item{iota_dist}{c("exponential", "weibull"); distribution of the infection period}
}
\value{
log likelihood
}
\description{
allows non-Markovian dynamics through Weibull-distribution infectious periods (Streftaris and Gibson, 2004)
allows generalized SIR process (Severo, 1969)
Severo, N. C. (1969). Generalizations of some stochastic epidemic models. Mathematical Biosciences, 4(3-4), 395-402.
Streftaris, G., & Gibson, G. J. (2004). Bayesian inference for stochastic epidemics in closed populations. Statistical Modelling, 4(1), 63-75.
}
