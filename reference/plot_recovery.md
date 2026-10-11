# Plot the recovery of true effects in simulated data

For simulated data, where the true effects are known, `plot_recovery()`
compares a fit with the truth. It shows the posterior intervals of the
sites together with the true effects, the proportion of intervals that
contain the true effect at each level, or the posterior means against
the true effects.

## Usage

``` r
plot_recovery(
  x,
  truth,
  type = c("intervals", "calibration", "scatter"),
  level = 0.9,
  n = 150,
  seed = NULL,
  scale = 1,
  ...
)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- truth:

  Numeric vector of true effects, one per site, in the same order as the
  sites in the data.

- type:

  Type of plot: `"intervals"` (the default), `"calibration"` or
  `"scatter"`. See Details.

- level:

  Posterior probability of the intervals, a number between 0 and 1. Only
  used with `type = "intervals"`.

- n:

  Maximum number of sites shown with `type = "intervals"`, an integer of
  at least 1.

- seed:

  Integer seed for the random choice of sites. If `NULL` (the default),
  the current state of R's random number generator is used.

- scale:

  Positive number by which the effects are multiplied in the plot, e.g.,
  `100` to show differences of proportions in percentage points. Not
  used with `type = "calibration"`.

- ...:

  These dots are for future extensions and must be empty.

## Value

A `ggplot` object.

## Details

`type` selects one of three plots.

- `"intervals"` shows the posterior interval and the true effect of each
  site, with the sites ordered by posterior mean, and marks the sites
  where the true effect lies outside the interval. If there are more
  than `n` sites, `n` of them are chosen at random; specify `seed` to
  make the choice reproducible.

- `"calibration"` shows, for levels from 5% to 95%, the proportion of
  all sites with the true effect inside the posterior interval of that
  level. The intervals are the ones that
  [`confint.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  returns. They are computed from the posterior distribution of each
  effect and not from a normal approximation, which would give wrong
  intervals for sites with a bimodal posterior distribution.

- `"scatter"` shows the posterior means against the true effects. The
  caption gives the root mean squared error of the posterior means and
  of the estimates.

The proportions in these plots describe one data set. They are not the
coverage of the intervals, which is the average of such a proportion
over repeated data sets. The intervals of all sites depend on the same
estimate of \\g\\, so they tend to miss together, and the proportion in
one data set varies around the nominal level more than it would for
independent sites.

## See also

[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
[`confint.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)

## Examples

``` r
truth <- twin_towers$theta_true
plot_recovery(twin_towers_fit, truth, seed = 1)

plot_recovery(twin_towers_fit, truth, type = "calibration")

plot_recovery(twin_towers_fit, truth, type = "scatter")

```
