% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{complete_theta}
\alias{complete_theta}
\title{Include beta to the parameter vector theta}
\usage{
complete_theta(theta, iota_dist, S0)
}
\arguments{
\item{theta}{parameters (without beta)}

\item{iota_dist}{distribution of infection periods}

\item{S0}{initial population size}
}
\value{
theta including beta
}
\description{
Include beta to the parameter vector theta
}
