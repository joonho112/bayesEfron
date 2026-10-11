# Summarize the posterior distributions of the site effects

`site_summary()` returns a data frame with one row per site: the site's
estimate and the mean, standard deviation, interval and mode of the
posterior distribution of its effect. Optionally, it also contains a
second, inner interval and the exceedance probability, i.e., the
posterior probability that the effect is greater than a threshold.
`exceedance_prob()` returns only that probability.

## Usage

``` r
site_summary(
  x,
  level = 0.9,
  inner_level = NULL,
  threshold = NULL,
  sites = NULL,
  scale = 1,
  reading = c("interpolated", "step")
)

exceedance_prob(
  x,
  threshold = 0,
  direction = c("above", "below"),
  sites = NULL,
  reading = c("interpolated", "step")
)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- level:

  Posterior probability of the interval in the columns `lower` and
  `upper`, a number between 0 and 1.

- inner_level:

  Optional level of a second interval, a number between 0 and 1. If
  specified, the columns `inner_lower` and `inner_upper` contain the
  limits of that interval. This is intended for plots that show two
  interval widths.

- threshold:

  Threshold for the exceedance probability, a number in the units of the
  estimates, whatever the value of `scale`. In `site_summary()` it is
  optional: if specified, the column `prob_above` contains the posterior
  probability that the effect is greater than `threshold`.

- sites:

  Sites to include, as a numeric vector of positions in the data or a
  character vector of site labels. If `NULL` (the default), all sites
  are included.

- scale:

  Positive number by which the effect columns are multiplied in the
  result, e.g., `100` to express differences in proportions as
  percentage points.

- reading:

  How quantiles and probabilities are computed for a fit with point
  masses: `"interpolated"` (the default) spreads the probability of each
  grid point evenly over one grid spacing around it, and `"step"` leaves
  it on the grid point. See the Details of
  [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md).

- direction:

  Either `"above"` (the default) for the probability that the effect is
  greater than `threshold`, or `"below"` for the probability that it is
  less than or equal to `threshold`.

## Value

`site_summary()` returns a data frame with one row per site and the
following columns.

- `site`, `label`:

  The site's position in the data and its label.

- `estimate`, `se`:

  The estimate and its standard error, as given.

- `mean`, `sd`:

  The posterior mean and the posterior standard deviation of the effect.

- `lower`, `upper`:

  The limits of the posterior interval at `level`.

- `mode`:

  The posterior mode.

- `sd_given_g`, `sd_about_g`:

  The two parts of `sd` described in Details.

- `inner_lower`, `inner_upper`:

  The limits of the interval at `inner_level`, only if `inner_level` is
  specified.

- `prob_above`:

  The posterior probability that the effect is greater than `threshold`,
  only if `threshold` is specified.

The values of `level`, `inner_level` and `threshold` are stored as
attributes.

`exceedance_prob()` returns a numeric vector of probabilities, named by
the site labels.

## Details

Both functions use the posterior distribution of each effect, which
averages over the posterior draws of \\g\\. They work for fits with or
without site-level draws (`site_draws` in
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)).

The posterior variance of an effect is the sum of two parts, and the
result contains the square root of each.

- `sd_given_g` is the square root of the variance that the effect would
  have if \\g\\ were known, averaged over the posterior draws of \\g\\.

- `sd_about_g` is the square root of the variance, over those draws, of
  the mean of the effect given \\g\\. This is the part that the
  uncertainty about \\g\\ adds.

The squares of the two columns sum to the square of `sd`. A plug-in
analysis, which treats one estimate of \\g\\ as known, omits the second
part. Its standard deviation is still not equal to `sd_given_g`, because
the posterior mean of \\g\\ is more spread out than a single draw of
\\g\\; for most sites it is larger than `sd_given_g`. See
[`as_normalmix()`](https://joonho112.github.io/bayesEfron/reference/as_normalmix.md)
for how to compute it.

For a fit with point masses, `reading` applies to the interval limits
and to the exceedance probabilities. With `reading = "step"`, a grid
point equal to `threshold` is counted as not above it, so the
probabilities for `direction = "above"` and `"below"` sum to 1.

The result does not include ranks: the rank of a site depends on the
effects of all sites jointly, and the effects are independent only given
\\g\\, so ranks cannot be computed from the separate posterior
distributions summarized here. Joint draws of the effects are needed for
that (see [`simulate()`](https://rdrr.io/r/stats/simulate.html) in
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)).

## See also

[`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md),
[`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md),
[`confint.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`summary.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)

## Examples

``` r
fit <- krw_fit

# The three firms with the largest posterior mean gaps, in percentage
# points, with the probability that the gap exceeds 2 points.
sites <- site_summary(fit, threshold = 0.02, scale = 100)
top <- sites[order(-sites$mean)[1:3], ]
top[, c("label", "estimate", "se", "mean", "lower", "upper", "prob_above")]
#>                        label estimate       se     mean    lower     upper
#> 62 Genuine Parts (Napa Auto) 9.817050 2.047973 6.414736 2.703855 11.605967
#> 28        Advance Auto Parts 7.396095 2.222589 3.862732 1.282151  7.228423
#> 61       O'Reilly Automotive 7.870053 2.441339 3.813371 1.136477  7.428633
#>    prob_above
#> 62  0.9840233
#> 28  0.8600060
#> 61  0.8370531

# How many firms have a positive gap with probability at least 0.95.
sum(exceedance_prob(fit) >= 0.95)
#> [1] 45

# The part of the posterior standard deviation that comes from not
# knowing g is largest for the firms with the largest estimates.
largest <- order(-sites$sd_about_g)[1:3]
sites[largest, c("label", "estimate", "se", "sd", "sd_given_g", "sd_about_g")]
#>                        label estimate       se       sd sd_given_g sd_about_g
#> 62 Genuine Parts (Napa Auto) 9.817050 2.047973 2.690078   2.275435  1.4348927
#> 61       O'Reilly Automotive 7.870053 2.441339 1.976101   1.847938  0.7000708
#> 28        Advance Auto Parts 7.396095 2.222589 1.849721   1.722103  0.6751528
```
