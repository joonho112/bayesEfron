# Package index

## Package

- [`bayesEfron-package`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  [`bayesEfron`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  : bayesEfron: Fully Bayesian Inference for the Empirical-Bayes
  Deconvolution Problem

## Fit

Primary entry point for fitting the fully Bayesian Efron log-spline
prior to a univariate random-effects deconvolution problem.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  : Fit the Bayesian Efron random-effects model

## Data preparation

Convert effect-size inputs (plain list or
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
object) into the canonical `bef_data` class consumed by the fit
pipeline.

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  : Convert input data to a bayesEfron data object

## Grid construction

Build the discrete support and natural-cubic-spline basis used by the
Stan model. Four recipes covering the paper’s real-data, simulation,
sensitivity, and an experimental KL-target setting.

- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  : Construct an Efron log-spline grid

## Compile and cache

Pre-warm the Stan model cache and maintain the on-disk and in-session
cache layers.

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  : Pre-compile the bayesEfron Stan model
- [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_clear_cache.md)
  : Clear bayesEfron compilation cache artifacts

## Diagnostics

Sampler-health and convergence diagnostics produced by the fit pipeline.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  : Extract sampler diagnostics from a bayesEfron fit

## Visualization

Caterpillar, prior-summary, sensitivity, and diagnostic views. Both
ggplot2 and base-graphics backends.

- [`plot(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  : Plot a bayesEfron random-effects fit

## S3 methods

Family-agnostic methods on the parent `bef_fit` class (`summary`,
`print`, `coef`, `vcov`, `confint`, `as.data.frame`, `nobs`, `logLik`,
[`posterior::as_draws`](https://mc-stan.org/posterior/reference/draws.html))
plus `format`/`print`/`summary` on `bef_data` and `bef_diagnostic`.

- [`format(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`format(`*`<summary.bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`format(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`format(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`print(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`summary(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`print(`*`<summary.bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`summary(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`coef(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`vcov(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`confint(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`as.data.frame(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`print(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`summary(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`print(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`summary(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`nobs(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`logLik(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  [`as_draws(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  : S3 methods for bayesEfron objects
