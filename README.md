
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PDSIR

<!-- badges: start -->
<!-- badges: end -->

The PDSIR R package implements an efficient data augmentation MCMC
(DA-MCMC) algorithm for fitting the general stochastic epidemic model to
incidence data. This package contains the code used in the paper
“Uniformly Ergodic Data-Augmented MCMC for Fitting the General
Stochastic Epidemic Model to Incidence Data” by R. Morsomme and J. Xu
available on ArXiv. The novelty of our DA-MCMC algorithm is the *joint*
update of the high-dimensional latent data in a Metropolis-Hastings
step. Compared to the existing single-site update DA-MCMC, our joint
proposal for the latent data significantly improves the mixing
properties of the resulting Markov chain.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rmorsomme/PDSIR")
```
