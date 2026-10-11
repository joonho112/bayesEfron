# Posterior predictive distribution for a new site

Draws from the posterior predictive distribution for a site that was not
in the data: either its true effect, or the estimate that a new study
with a given standard error would produce. For a new site whose estimate
is already known, [`predict()`](https://rdrr.io/r/stats/predict.html)
instead returns the mean, standard deviation, interval and mode of the
posterior distribution of its effect.

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
  seed = NULL,
  level = 0.9
)
```

## Arguments

- object:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- newdata:

  For `type = "theta_hat"`, a numeric vector of positive standard
  errors, one per new site. It is required for that type and not used
  for `type = "theta"`. Can also be a data frame with the columns
  `estimate` and `se`, for new sites whose estimates are known; see
  Details.

- ...:

  These dots are for future extensions and must be empty.

- type:

  Quantity to draw for a new site: `"theta"` (the default) for its true
  effect, or `"theta_hat"` for its estimate.

- n:

  Number of draws, a positive integer. If `NULL` (the default), the
  number of posterior draws in the fit is used. If `n` is larger than
  that, posterior draws of \\g\\ are used more than once.

- jitter:

  If `TRUE`, the draws from a fit with point masses are spread between
  the grid points; see Details.

- seed:

  Integer seed for the draws. If `NULL` (the default), no seed is set.
  If a seed is specified, the state of R's random number generator is
  restored afterward.

- level:

  Posterior probability of the interval in the columns `lower` and
  `upper`, a number between 0 and 1. Only used if `newdata` is a data
  frame.

## Value

For `type = "theta"`, and for `type = "theta_hat"` with one standard
error in `newdata`, a numeric vector of `n` draws. For
`type = "theta_hat"` with more than one standard error, a matrix with
`n` rows and one column per standard error.

If `newdata` is a data frame, a data frame with one row per row of
`newdata` and the following columns.

- `estimate`, `se`:

  The estimate and its standard error, as given.

- `mean`, `sd`:

  The posterior mean and the posterior standard deviation of the effect.

- `lower`, `upper`:

  The limits of the posterior interval at `level`.

- `mode`:

  The posterior mode.

## Details

[`predict()`](https://rdrr.io/r/stats/predict.html) has three uses,
which are selected by `type` and `newdata`.

- `type = "theta"` (the default) draws the true effect of a new site
  from \\g\\. Each draw uses a different posterior draw of \\g\\, so the
  draws describe the distribution of effects and include the uncertainty
  about \\g\\.

- `type = "theta_hat"` draws the estimate that a new study would
  produce: sampling error with the standard error given in `newdata` is
  added to a draw of the effect. If `newdata` contains several standard
  errors, the result has one column per new site. Each row then uses one
  posterior draw of \\g\\ for all sites and draws a separate effect and
  sampling error for each site, so the sites in a row are independent
  given that draw.

- If `newdata` is a data frame with the columns `estimate` and `se`, the
  estimates of the new sites are known. Nothing is drawn, and for each
  row [`predict()`](https://rdrr.io/r/stats/predict.html) returns the
  posterior mean, the posterior standard deviation, the interval at
  `level` and the posterior mode of the effect. `type`, `n`, `jitter`
  and `seed` are not used.

For a fit with point masses, \\g\\ is a distribution on the grid, so
each draw of an effect is one of the grid points. If `jitter = TRUE`, a
uniform random number between minus and plus half the grid spacing is
added to the draw, which spreads the draws over the intervals between
the grid points. For a fit with `component = "normal"`, a draw comes
from the normal component at the grid point that was drawn, and `jitter`
is not used.

A new site with a known estimate was not in the data, so before its
estimate is observed the distribution of its effect is the posterior
mean of \\g\\. Bayes' rule with the estimate then gives the posterior
distribution that [`predict()`](https://rdrr.io/r/stats/predict.html)
summarizes. The summaries are those that the site would have in a fit of
the data with that site added, on the same grid, up to Monte Carlo
error. They differ from the summaries that
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
gives for a site in the fit with the same estimate and standard error,
because that site's estimate was already used to estimate \\g\\. The
difference is small when there are many sites and is largest for an
estimate in a tail of the distribution. Each row of `newdata` is treated
as the only new site, so the other rows do not affect its result. For a
fit with point masses, the interval limits are computed as in
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
with its default `reading`.

## See also

[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
# The callback gap of a firm that was not in the data.
theta_new <- predict(krw_fit, seed = 1)
quantile(theta_new, c(0.05, 0.5, 0.95))
#>           5%          50%          95% 
#> -0.008055305  0.013672701  0.045057598 
mean(theta_new > 0)
#> [1] 0.85625

# An estimate with a standard error of 1 percentage point.
estimate_new <- predict(krw_fit, newdata = 0.01, type = "theta_hat",
                        seed = 1)
quantile(estimate_new, c(0.05, 0.5, 0.95))
#>          5%         50%         95% 
#> -0.01264502  0.01529825  0.04897299 

# Two new firms whose estimates are known: gaps of 0 and of 5 percentage
# points, with standard errors of 1 and 2 points.
predict(krw_fit, newdata = data.frame(estimate = c(0, 0.05),
                                      se = c(0.01, 0.02)))
#>   estimate   se        mean          sd        lower      upper        mode
#> 1     0.00 0.01 0.005452166 0.007714249 -0.007174759 0.01826861 0.006430032
#> 2     0.05 0.02 0.028766853 0.014235719  0.007674624 0.05445411 0.025743815
```
