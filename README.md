
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SSMimpute

<!-- badges: start -->
<!-- badges: end -->

The goal of SSMimpute is to …

## Installation

You can install the development version of SSMimpute like so:

``` r
devtools::install_github("js5959/SSMimpute")
```

The workflow is: - Data exploration: Using
``` shrinkTVP()``explore_SSM ``` function to get some sense of the data,
identify the formula, check the possible ARIMA model of latent space and
find the change points. The problem of over fitting is partially
alleviated by introducing shrinkage. Based on the exploration result, we
have the
