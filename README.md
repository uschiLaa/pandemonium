
# pandemonium

<!-- badges: start -->
<!-- badges: end -->

The goal of pandemonium is to interactively explore clustering solutions
of physics problems in both observable and parameter space in a Shiny app.

## Installation

You can install the package from GitHub using devtools:

``` r
# install.packages("devtools")
devtools::install_github("uschiLaa/pandemonium")
```

## Example

The package comes with an example data set `b_anomaly` which contains the three
necessary inputs to the app: theory predictions and experimentally observed
values for all observables (pred and exp), full inverse covariance matrix
(covInv) and parameter values (wc).

Using the example data the app is launched as:

``` r
library(pandemonium)
pandemonium(b_anomaly$pred, b_anomaly$covInv, b_anomaly$wc, b_anomaly$exp)
```
