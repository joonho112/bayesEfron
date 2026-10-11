# Plot a fitted model

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
a fitted model shows the plot selected by `type`. Each type is produced
by one of the package's plot functions, e.g., `"sites"` by
[`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md)
and `"g"` by
[`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
and further arguments are passed on to that function. See Details for
all types.

## Usage

``` r
# S3 method for class 'bef_fit_re'
plot(x, type = "sites", ..., backend = NULL)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- type:

  The plot to show, e.g., `"sites"` (the default), `"g"` or
  `"shrinkage"`. See Details for all types.

- ...:

  Further arguments passed to the plot function for `type`, such as
  `level`, `sort_by` or `scale`.

- backend:

  Deprecated. `backend = "base"` plots the types `"caterpillar"`, `"g"`,
  `"diagnostic"` and `"loo"` with base graphics instead of ggplot2 and
  issues a warning. This option will be removed.

## Value

A `ggplot` object, which can be modified with ggplot2 functions, e.g.,
by adding layers. For `type = "mcmc"`, the value returned by
[`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md).
With `backend = "base"`, the plot is shown on the current graphics
device and `NULL` is returned invisibly.

## Details

`type` can be one of the following. The help page of each function
describes the plot and its arguments.

- `"sites"` (the default):
  [`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md),
  a caterpillar plot of the site effects with their posterior intervals

- `"g"`:
  [`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
  the estimated distribution of effects

- `"shrinkage"`:
  [`plot_shrinkage()`](https://joonho112.github.io/bayesEfron/reference/plot_shrinkage.md),
  a shrinkage plot of posterior means against estimates

- `"exceedance"`:
  [`plot_exceedance()`](https://joonho112.github.io/bayesEfron/reference/plot_exceedance.md),
  the probability that each effect exceeds a threshold

- `"site_posterior"`:
  [`plot_site_posterior()`](https://joonho112.github.io/bayesEfron/reference/plot_site_posterior.md),
  the posterior distribution of the effect of each site specified in
  `sites`

- `"diagnostics"`:
  [`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md),
  sampler and grid diagnostics

- `"loo"`:
  [`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md),
  the Pareto k values of the sites from leave-one-out cross-validation

- `"ppc"`:
  [`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md),
  a posterior predictive check of the model against the estimates

- `"mcmc"`:
  [`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md),
  trace plots and other plots of the posterior draws

- `"grid"`:
  [`plot.bef_grid()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid.md),
  the grid in relation to the estimates

- `"recovery"`:
  [`plot_recovery()`](https://joonho112.github.io/bayesEfron/reference/plot_recovery.md),
  a comparison of the fit with the known true effects specified in
  `truth`

`"caterpillar"` and `"diagnostic"` are accepted as aliases for `"sites"`
and `"diagnostics"`.

`type = "g"` shows the density of \\g\\. To plot the probabilities of
the grid points instead, call `plot_g(x, type = "mass")`. Note that
their heights depend on the number of grid points.

## See also

[`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md),
[`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`vignette("plots")`](https://joonho112.github.io/bayesEfron/articles/plots.md)

## Examples

``` r
plot(krw_fit, n_extreme = 10, scale = 100)

plot(krw_fit, type = "g", scale = 100)

plot(krw_fit, type = "diagnostics")


# The result is a ggplot object, so labels and layers can be added.
plot(krw_fit, type = "shrinkage", scale = 100) +
  ggplot2::labs(x = "Estimated gap (percentage points)")

```
