# Trace plots and other MCMC plots with bayesplot

`plot_mcmc()` is a convenience wrapper around bayesplot functions for
the posterior draws of a fit: trace plots, rank plots and densities by
chain, pairs plots, and the energy plot of the sampler. For any other
plot, pass the draws (`x$draws`) and `bayesplot::nuts_params(x)` to
bayesplot directly.

## Usage

``` r
plot_mcmc(
  x,
  type = c("trace", "rank", "dens", "pairs", "energy"),
  pars = c("alpha", "lambda", "mean_g", "sd_g"),
  ...
)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- type:

  Type of plot: `"trace"` (the default), `"rank"`, `"dens"`, `"pairs"`
  or `"energy"`. See Details.

- pars:

  Character vector with the names of the parameters to plot. A name
  without an index selects all elements of that parameter, e.g.,
  `"alpha"` selects every spline coefficient. By default the spline
  coefficients, `lambda`, and the mean and standard deviation of \\g\\
  are plotted. Not used with `type = "energy"`.

- ...:

  These dots are for future extensions and must be empty.

## Value

The object returned by the bayesplot function: a `ggplot` object or, for
`type = "pairs"`, a grid of plots that is shown when it is printed.

## Details

Each type calls one bayesplot function: `"trace"` calls `mcmc_trace()`,
`"rank"` calls `mcmc_rank_overlay()`, `"dens"` calls
`mcmc_dens_overlay()`, `"pairs"` calls `mcmc_pairs()` (which also marks
the divergent transitions), and `"energy"` calls `mcmc_nuts_energy()`.
The bayesplot package must be installed, and the plots use its theme and
colors. `type = "pairs"` also needs the gridExtra package.

Only the draws stored in the fit are plotted. If the fit was thinned
(`thin`), these are a subset of the iterations, whereas
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
and
[`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md)
describe the whole run. The energy plot compares the energy at each
iteration with its change from the previous iteration, so it needs every
iteration. For a thinned fit, including the fits stored in the package,
`type = "energy"` therefore throws an error. In that case
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
still reports the E-BFMI of the whole run.

In a trace plot, chains that have converged cover the same range and
none of them drifts. See
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
and
[`vignette("plots")`](https://joonho112.github.io/bayesEfron/articles/plots.md)
for how to read the trace plot and the rank plot.

## See also

[`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)

## Examples

``` r
plot_mcmc(krw_fit, pars = c("lambda", "mean_g", "sd_g"))

plot_mcmc(krw_fit, type = "rank", pars = "alpha")

plot_mcmc(krw_fit, type = "dens", pars = c("mean_g", "sd_g"))
```
