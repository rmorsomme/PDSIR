
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PDSIR

<!-- badges: start -->
<!-- badges: end -->

The `PDSIR` R package implements an efficient data augmentation MCMC
(DA-MCMC) algorithm for exact Bayesian inference under the semi-Markov
stochastic susceptible-infectious-removed (SIR) model, given discretely
observed counts of infections. The novelty of this DA-MCMC algorithm is
the *joint* update of the high-dimensional latent data. In a
Metropolis-Hastings step, the latent data are jointly proposed from a
surrogate process carefully designed to closely resemble the target
process and from which we can efficiently generate epidemics consistent
with the observed data. This yields a MCMC algorithm that explores the
high-dimensional latent space efficiently, mixes significantly better
than single-site samplers, and scales to outbreaks with thousands of
infections.

## Installation

You can install the development version of PDSIR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rmorsomme/PDSIR")
```

## Illustration of synthetic data

We employ the DA-MCMC on artifical data; semi-Marko model

This package contains the code used in the paper “Uniformly Ergodic
Data-Augmented MCMC for Fitting the General Stochastic Epidemic Model to
Incidence Data” by R. Morsomme and J. Xu available on ArXiv. We use it
to fit a semi-Markov susceptible-infectious-removed model to the
2013-2015 outbreak of Ebola Haemorrhagic Fever in Gu'eck'edou, Guinea.

This is a basic example which shows you how to solve a common problem:

``` r
library(PDSIR)

# Setup
set.seed(0)
S0 <- 500  # initial number of susceptible individuals
I0 <- 5    # initial number of infectious individuals

t_end <- 6 # # duration of observation period

iota_dist <- "weibull" # distribution of the infection periods
theta <- list(
  R0 = 2,               # basic reproduction number
  lambda = 1, shape = 2 # parameters of the Weibull distribution for the infection periods
  )  

theta <- complete_theta(theta, iota_dist, S0) # add the infection rate parameter, beta, and the average infection period, gamma

# Simulate artificial data
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

# Trajectories of compartments
draw_trajectories(SIR, t_end)
```

<img src="man/figures/README-synthetic-data-1.png" width="100%" />

Observed data

``` r
K <- 10 # number of intervals
Y <- observed_data(SIR, K)

print(Y$ts ) # endpoints of intervals
#>  [1] 0.0 0.6 1.2 1.8 2.4 3.0 3.6 4.2 4.8 5.4 6.0
print(Y$T_k) # number of infections per interval
#>  [1]  7 11 38 51 67 88 63 39 18  8
```

## Illustration on the 2013-2015 outbreak of Ebola in Western Africa

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
