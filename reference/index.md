# Package index

## Fit the model

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
fits the model, and
[`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
describes the object that it returns and the methods for it.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  : Fit the Efron log-spline model to estimates and standard errors
- [`print(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`summary(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`summary(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`print(`*`<summary.bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`coef(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`confint(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`vcov(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`fitted(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`as.data.frame(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`nobs(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`logLik(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`as_draws(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`format(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`format(`*`<summary.bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  : Fitted Efron log-spline model objects
- [`predict(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md)
  : Predict the effect or the estimate at a new site

## Prepare the data and the grid

Check the estimates and their standard errors, and build the grid and
the spline basis on which the distribution of effects is represented.

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`print(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`summary(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`format(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  : Extract estimates and standard errors for a fit
- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  : Build the grid and spline basis for the effect distribution

## Check and plot a fit

Convergence diagnostics and plots of a fitted model.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`print(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`summary(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`format(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  : Convergence diagnostics of a fitted model
- [`plot(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  : Plot a fitted model

## Compare models

Leave-one-out cross-validation for one fit, and a comparison of several
specifications of the grid and the spline basis.

- [`loo(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  [`waic(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  : Leave-one-out cross-validation for a fitted model
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  [`print(`*`<bef_grid_comparison>`*`)`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  : Compare grids and numbers of spline functions by cross-validation

## Compile the model

Compile the Stan program before the first fit, or delete the compiled
model.

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  : Compile the Stan model ahead of time, or remove the compiled model

## Data

The teacher expectancy studies and the model fitted to them.

- [`raudenbush1985`](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
  : Teacher expectancy studies
- [`raudenbush_fit`](https://joonho112.github.io/bayesEfron/reference/raudenbush_fit.md)
  : Model fitted to the teacher expectancy studies

## Package

- [`bayesEfron-package`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  [`bayesEfron`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  : bayesEfron: Fully Bayesian Inference for the Empirical-Bayes
  Deconvolution Problem
