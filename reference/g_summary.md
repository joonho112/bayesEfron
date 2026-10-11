# Summarize the estimated distribution of effects

Computes the posterior mean and a pointwise posterior interval for a
feature of \\g\\, the distribution of true effects: its density, its
distribution function, its probabilities on the grid, its quantiles, the
probability of an interval of effects, or its mean, variance and
standard deviation. The feature is computed for each posterior draw of
\\g\\ and then summarized over the draws.

## Usage

``` r
g_summary(
  x,
  type = c("density", "cdf", "mass", "quantile", "prob", "moments"),
  level = 0.9,
  at = NULL,
  probs = NULL,
  lower = -Inf,
  upper = Inf,
  scale = 1,
  reading = c("interpolated", "step")
)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- type:

  Feature of \\g\\ to summarize: `"density"` (the default), `"cdf"` for
  the distribution function, `"mass"` for the probabilities at the grid
  points, `"quantile"`, `"prob"` for the probability of an interval of
  effects, or `"moments"` for the mean, variance and standard deviation.
  The standard deviation of \\g\\ measures the heterogeneity of the
  effects.

- level:

  Posterior probability of the pointwise interval, a number between 0
  and 1.

- at:

  Numeric vector of effect values at which the density or the
  distribution function is evaluated. If `NULL` (the default), the grid
  points are used.

- probs:

  Numeric vector of probabilities for `type = "quantile"`, each greater
  than 0 and less than 1. If `NULL` (the default), 0.05, 0.25, 0.5, 0.75
  and 0.95 are used.

- lower, upper:

  Limits of the interval of effects for `type = "prob"`. The result is
  the probability that an effect is greater than `lower` and less than
  or equal to `upper`.

- scale:

  Positive number by which effect values are multiplied in the result,
  e.g., `100` to express differences in proportions as percentage
  points. See Details.

- reading:

  How quantiles and probabilities are computed for a fit with point
  masses. With `"interpolated"` (the default), the probability of each
  grid point is spread evenly over one grid spacing around it. With
  `"step"`, it stays on the grid point. Use `"step"` if the effects must
  stay within given bounds that were specified as the limits of the grid
  (see `support` in
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)).
  See Details.

## Value

A data frame with the following columns. The interval level is stored in
the attribute `level`.

- `value`:

  The effect value, for `"density"`, `"cdf"` and `"mass"`.

- `prob`:

  The probability, for `"quantile"`.

- `from`, `to`:

  The limits of the interval of effects, for `"prob"`.

- `parameter`:

  The name of the moment, for `"moments"`.

- `mean`:

  The posterior mean.

- `lower`, `upper`:

  The limits of the pointwise posterior interval.

## Details

The model estimates \\g\\ as probabilities on the grid points. For a fit
with point masses (`component = "point"`, the default), there are two
ways to compute a quantile or a probability from them, and `reading`
selects one.

- `"interpolated"` spreads the probability of each grid point evenly
  over an interval of one grid spacing centered at the point. The
  distribution function is then continuous, and the results change
  little when the number of grid points changes. For the first and the
  last grid point this interval extends half a grid spacing beyond the
  limits of the grid, so a quantile can lie up to half a spacing outside
  the grid, and a small amount of probability lies beyond each limit.

- `"step"` leaves the probability on the grid point, as the model
  specifies, so all results stay on the grid. The distribution function
  then jumps at the grid points, and a probability such as \\P(\theta \>
  0)\\ depends on where the grid points lie relative to zero.

For `type = "density"`, the probabilities of point masses are always
divided by the grid spacing, which corresponds to `"interpolated"`.

The posterior interval is pointwise: at each value it has posterior
probability `level` for that value, not for the whole curve at once.

`at`, `lower` and `upper` are specified in the units of the estimates.
`scale` only changes the units of the result: effect values are
multiplied by it, densities are divided by it, and probabilities are
unchanged. For example, if the estimates are differences in proportions,
`scale = 100` gives percentage points.

## See also

[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
[`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
[`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
[`confint.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)

## Examples

``` r
fit <- krw_fit

# The proportion of firms whose callback gap favors white applicants,
# and the proportion with a gap above 2 percentage points.
g_summary(fit, type = "prob", lower = 0)
#>   from  to      mean     lower     upper
#> 1    0 Inf 0.8754888 0.8083673 0.9366439
g_summary(fit, type = "prob", lower = 0.02)
#>   from  to      mean     lower     upper
#> 1 0.02 Inf 0.3641116 0.2582442 0.4805717

# Quantiles of the gaps, in percentage points.
g_summary(fit, type = "quantile", probs = c(0.05, 0.5, 0.95), scale = 100)
#>   prob       mean     lower      upper
#> 1 0.05 -0.6379926 -1.275347 -0.1101583
#> 2 0.50  1.5080321  1.131837  1.9253976
#> 3 0.95  4.4982003  3.385929  6.2305993

# The density at the grid points, which is what plot_g() shows.
head(g_summary(fit, scale = 100))
#>       value         mean        lower       upper
#> 1 -8.289621 0.0018215418 2.027327e-06 0.007353329
#> 2 -8.048199 0.0012741173 1.945940e-06 0.005158977
#> 3 -7.806777 0.0009383640 1.526605e-06 0.003873485
#> 4 -7.565355 0.0007200636 1.715523e-06 0.002983518
#> 5 -7.323932 0.0005729214 1.246803e-06 0.002355854
#> 6 -7.082510 0.0004715501 9.368827e-07 0.001985074
```
