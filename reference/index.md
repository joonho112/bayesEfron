# Package index

## Fit the model

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
fits the model.
[`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
describes the fitted object and its methods.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  : Fit the Bayesian Efron log-spline model to estimates and standard
  errors
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
  [`quantile(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`format(`*`<summary.bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`residuals(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`simulate(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`ndraws(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`variables(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`AIC(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  [`BIC(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  : Fitted Efron log-spline model objects
- [`predict(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md)
  : Posterior predictive distribution for a new site

## Prepare the data and the grid

Check the estimates and their standard errors, and build the grid and
the spline basis on which the distribution of effects is represented.

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`print(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`summary(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`format(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  [`as.data.frame(`*`<bef_data>`*`)`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  : Extract estimates and standard errors from a data frame or list
- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  [`format(`*`<bef_grid>`*`)`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  [`print(`*`<bef_grid>`*`)`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  : Build the grid and spline basis for the distribution of effects

## Summarize a fit

The estimated distribution of effects and the posterior distributions of
the site effects, as data frames.

- [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md)
  : Summarize the estimated distribution of effects
- [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  [`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  : Summarize the posterior distributions of the site effects

## Check a fit

Sampler and grid diagnostics.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`print(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`summary(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`as.data.frame(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`format(`*`<bef_diagnostic>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`nuts_params(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`log_posterior(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`rhat(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  [`neff_ratio(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  : MCMC convergence and grid diagnostics for a fitted model

## Plots

Each plot function returns a `ggplot` object.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit calls
them by `type`.

- [`plot(`*`<bef_fit_re>`*`)`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  : Plot a fitted model

### Results

- [`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md)
  : Density and CDF plots of the estimated distribution of effects
- [`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md)
  : Caterpillar (forest) plot of the site effects
- [`plot_shrinkage()`](https://joonho112.github.io/bayesEfron/reference/plot_shrinkage.md)
  : Shrinkage plot of posterior means against estimates
- [`plot_exceedance()`](https://joonho112.github.io/bayesEfron/reference/plot_exceedance.md)
  : Plot the exceedance probabilities of the site effects
- [`plot_site_posterior()`](https://joonho112.github.io/bayesEfron/reference/plot_site_posterior.md)
  : Plot the posterior distribution of individual site effects

### Checks of a fit

- [`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md)
  : Plot sampler diagnostics and grid diagnostics
- [`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md)
  : Plot the Pareto k values of leave-one-out cross-validation
- [`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md)
  : Posterior predictive check
- [`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md)
  : Trace plots and other MCMC plots with bayesplot

### Grids and simulated data

- [`plot(`*`<bef_grid>`*`)`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid.md)
  : Plot a grid and its spline basis with the estimates
- [`plot(`*`<bef_grid_comparison>`*`)`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid_comparison.md)
  : Plot the sensitivity of the results to the grid
- [`plot_recovery()`](https://joonho112.github.io/bayesEfron/reference/plot_recovery.md)
  : Plot the recovery of true effects in simulated data

## Compare models

Leave-one-out cross-validation, a comparison of fits with different
grids and spline bases, and the estimate of g in the form that the
empirical Bayes functions of the ebnm package take.

- [`loo(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  [`waic(`*`<bef_fit>`*`)`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  : Leave-one-out cross-validation and WAIC for a fitted model
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  [`print(`*`<bef_grid_comparison>`*`)`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  : Compare fits on several grids (sensitivity analysis)
- [`as_normalmix()`](https://joonho112.github.io/bayesEfron/reference/as_normalmix.md)
  : Export the estimated distribution of effects for the ebnm package

## Compile the model

Compile the Stan program before the first fit, or delete the compiled
model.

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  : Compile the Stan programs before the first fit

## Data sets and stored fits

Example data, and fits of them that are stored in the package so that
the examples run without CmdStan.

### Callback gaps of 97 firms, and 1,500 simulated sites (Lee and Sui, 2025)

- [`krw_discrimination`](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md)
  : Racial callback gaps among 97 large firms
- [`krw_fit`](https://joonho112.github.io/bayesEfron/reference/krw_fit.md)
  : Model fitted to the racial callback gaps
- [`twin_towers`](https://joonho112.github.io/bayesEfron/reference/twin_towers.md)
  : Twin Towers simulated data
- [`twin_towers_fit`](https://joonho112.github.io/bayesEfron/reference/twin_towers_fit.md)
  : Model fitted to the Twin Towers simulated data

### A small meta-analysis (19 teacher expectancy studies)

- [`raudenbush1985`](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
  : Teacher expectancy studies
- [`raudenbush_fit`](https://joonho112.github.io/bayesEfron/reference/raudenbush_fit.md)
  : Model fitted to the teacher expectancy studies

### A multisite trial in education (79 schools of Project STAR)

- [`star_schools`](https://joonho112.github.io/bayesEfron/reference/star_schools.md)
  : Small-class effects in 79 schools of Project STAR

## Package

- [`bayesEfron-package`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  [`bayesEfron`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-package.md)
  : bayesEfron: Fully Bayesian Deconvolution with Efron's Log-Spline
  Prior
