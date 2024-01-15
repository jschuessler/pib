# pib - Sensitivity analysis for post-instrument bias

This R package implements the tools described in:

Schuessler, J., Glynn, A. N., & Rueda, M.R. (2023). "Post-Instrument Bias".

<!-- Glynn, A. N., Rueda, M. R., & Schuessler, J. (2023). Post-Instrument Bias in Linear Models. Sociological Methods & Research. -->

## Installation

``` r
if(!require(devtools)) install.packages("devtools") 
library(devtools) 
devtools::install_github("jschuessler/pib")
```

## Usage

The workflow consists of three steps:

1.  `estimate_sensitivity` takes as input a dataframe and the name of the outcome, treatment, instrument, post-instrument, and control variables. The function computes bootstrap replicates for the required model components for the sensitivity analysis.
2.  `post_instrument_bounds` computes bounds given the output from `estimate_sensitivity` and values for the sensitivity parameters `gamma` and (optionally) `sigma`
3.  `plot_pib_bounds` is used to plot the output from `post_instrument_bounds`

``` r
out <- estimate_sensitivity(df, outcome = y, treat = d, inst = z, post.inst = m)

gamma <- c(-0.02, 0.02)
sigma <- c(0, 0.002)

bounds <- post_instrument_bounds(gamma = gamma, sigma = sigma,
                                 estimates = out)

plot_pib_bounds(bounds)
```

For an extended example, see [here](https://github.com/jschuessler/pib/blob/main/vignettes/pib-replication.Rmd).
