% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ebola_gueckedou}
\alias{ebola_gueckedou}
\title{Weekly infection counts in Gueckedou}
\format{
## `ebola_gueckedou`
A list of three vectors:
\describe{
  \item{I_k}{number of weekly infections}
  \item{dates}{date of the week}
  \item{ts}{endpoints of observation intervals (units: days), used in the MCMC algorithm; the intervals correspond to the weeks in `dates`}
}
}
\usage{
ebola_gueckedou
}
\description{
Weekly infection counts of Ebola in the prefecture Gueckedou in Guinea during the 2013-2025 outbreak.
The full data for Guinea can be found in `ebola_guinea`.
}
\keyword{datasets}
