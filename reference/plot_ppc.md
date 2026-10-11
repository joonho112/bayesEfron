# Posterior predictive check

`plot_ppc()` is a graphical posterior predictive check. It compares the
estimates with the distribution that the fitted model predicts for them,
through the probability integral transform (PIT) of each estimate or
through the marginal density of the estimates.

## Usage

``` r
plot_ppc(
  x,
  type = c("pit_ecdf", "pit_hist", "marginal"),
  loo = FALSE,
  level = 0.95,
  scale = 1,
  ...
)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- type:

  Type of plot: `"pit_ecdf"` (the default), `"pit_hist"` or
  `"marginal"`. See Details.

- loo:

  If `TRUE`, the PIT of each site is computed from the posterior
  distribution without that site (leave-one-out). Not used with
  `type = "marginal"`.

- level:

  Level of the pointwise reference intervals, a number between 0 and 1.
  Not used with `type = "marginal"`.

- scale:

  Positive number by which the estimates are multiplied on the
  horizontal axis. Only used with `type = "marginal"`.

- ...:

  These dots are for future extensions and must be empty.

## Value

A `ggplot` object.

## Details

Under the model, the estimate of a site is an effect from \\g\\ plus
normal sampling error, so its distribution is \\g\\ convolved with a
normal distribution whose standard deviation is the site's standard
error. The PIT of an estimate is the probability that the model assigns
to a smaller estimate at that site. If \\g\\ were known and the model
were correct, the PIT values would be uniformly distributed between 0
and 1, so the check compares the PIT values of the fitted model with the
uniform distribution.

For site \\i\\ with estimate \\\hat\theta_i\\ and standard error
\\\sigma_i\\, and for one posterior draw of the probabilities \\g_j\\ at
the grid points \\\tau_j\\, \$\$\mathrm{PIT}\_i = \sum_j g_j \\
\Phi\left(\frac{\hat\theta_i - \tau_j} {\sqrt{\sigma_i^2 +
\tau_c^2}}\right),\$\$ where \\\Phi\\ is the standard normal
distribution function and \\\tau_c\\ is the standard deviation of a
normal component (zero for point masses). The plots show the mean of
\\\mathrm{PIT}\_i\\ over the posterior draws.

There are three types of plot.

- `"pit_ecdf"` shows the empirical distribution function of the PIT
  values minus the distribution function of the uniform distribution,
  with pointwise reference intervals for independent uniform values.

- `"pit_hist"` shows a histogram of the PIT values, with the
  corresponding reference intervals for the height of a bar. If there
  are fewer than 30 sites, the bars would contain too few values, and
  the `"pit_ecdf"` plot is shown instead.

- `"marginal"` shows a histogram of the estimates together with the
  density that the model assigns to an estimate, averaged over the
  standard errors of the sites.

The reference intervals are computed for independent uniform values,
which is what the PIT values would be if \\g\\ were known. However, the
PIT values of a fitted model all depend on the same estimate of \\g\\,
and the band ignores this dependence. Its level is therefore not the
probability that the line of a fitted model lies inside the band at a
given point, with or without `loo = TRUE`.

Because the model was fitted to the same estimates that it is checked
against, the check is lenient. With `loo = TRUE`, the PIT of each site
is computed from the posterior distribution without that site, using
Pareto smoothed importance sampling. This requires the loo package.

The check can show that the estimates are not what a normal likelihood
with the given standard errors and some \\g\\ on this grid would
produce, e.g., because the standard errors are too small, the grid is
too narrow, or there are too few spline functions. Note that it cannot
show that the shape of \\g\\ is correct. Distributions of effects with
different shapes give nearly the same distribution of noisy estimates,
so a model can pass the check with a wrong \\g\\.

## See also

[`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md),
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md),
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)

## Examples

``` r
plot_ppc(krw_fit)

plot_ppc(twin_towers_fit, type = "pit_hist")

plot_ppc(krw_fit, type = "marginal", scale = 100)

```
