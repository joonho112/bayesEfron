# Predict the effect or the estimate at a new site

[`predict()`](https://rdrr.io/r/stats/predict.html) draws from the
predictive distribution for a site that was not in the data. With
`type = "theta"` it draws the true effect of the new site from \\g\\,
using a different posterior draw of \\g\\ each time, so that the draws
describe the distribution of effects with the uncertainty about \\g\\
included. With `type = "theta_hat"` it adds sampling error with the
standard error given in `newdata`, and so draws the estimate that a new
study of that precision would give.

## Usage

``` r
# S3 method for class 'bef_fit_re'
predict(
  object,
  newdata = NULL,
  ...,
  type = c("theta", "theta_hat"),
  n = NULL,
  jitter = FALSE,
  seed = NULL
)
```

## Arguments

- object:

  A `bef_fit_re` object from
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- newdata:

  Numeric vector of positive standard errors, one for each new site. It
  is needed for `type = "theta_hat"` and not used otherwise.

- ...:

  These dots are for future extensions and must be empty.

- type:

  Character string: `"theta"` (the default) for the true effect of a new
  site, or `"theta_hat"` for its estimate.

- n:

  Whole number, the number of draws. `NULL` (the default) means the
  number of posterior draws in the fit. If `n` is larger than that,
  posterior draws of \\g\\ are used more than once.

- jitter:

  Logical. If `TRUE`, the draws are spread between the grid points, as
  described in Details. Defaults to `FALSE`.

- seed:

  Whole number, a seed for the draws, or `NULL` (the default). When a
  seed is given, R's random number generator is left afterward in the
  state in which it was found.

## Value

A numeric vector of `n` draws. For `type = "theta_hat"` with more than
one standard error in `newdata`, a matrix with `n` rows and one column
for each standard error.

## Details

\\g\\ is a distribution on the grid, so a draw of the effect is one of
the grid points. With `jitter = TRUE` a uniform random number of up to
half the grid spacing is added or subtracted, which spreads the draws
over the intervals between the grid points.

## See also

[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
# The effect in a new study.
theta_new <- predict(raudenbush_fit, seed = 1)
quantile(theta_new, c(0.05, 0.5, 0.95))
#>    5%   50%   95% 
#> -0.23  0.07  0.73 
mean(theta_new > 0)
#> [1] 0.67625

# The estimate that a new study with standard error 0.15 would give.
estimate_new <- predict(raudenbush_fit, newdata = 0.15, type = "theta_hat",
                        seed = 1)
quantile(estimate_new, c(0.05, 0.5, 0.95))
#>          5%         50%         95% 
#> -0.33594410  0.07391744  0.76314394 
```
