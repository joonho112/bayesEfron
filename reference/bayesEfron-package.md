# bayesEfron: Fully Bayesian Deconvolution with Efron's Log-Spline Prior

Estimates the distribution of true effects from a set of estimates with
known standard errors, as in a meta-analysis or a multisite trial, using
the fully Bayesian version of Efron's log-spline prior described in Lee
and Sui (2025)
[doi:10.3390/math13162639](https://doi.org/10.3390/math13162639) . The
fit returns posterior draws of the distribution of effects, from which
the summaries of each effect are computed. The package also provides
sampler and grid diagnostics, plots made with 'ggplot2', and model
comparison by approximate leave-one-out cross-validation. The estimated
distribution can be passed to the 'ebnm' package for comparison with
empirical Bayes estimates. The model is written in 'Stan' and is fitted
with 'CmdStan' through the 'cmdstanr' package.

## Details

To get started, read
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md),
which walks through a complete analysis of the callback gaps of 97
firms. The model, its prior distributions and the way it is fitted are
described in
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

## Functions

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  extracts the estimates and their standard errors from a list, a data
  frame or a
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  object and checks them.

- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  builds the grid and the spline basis, and
  [plot()](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid.md)
  shows them together with the estimates.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  fits the model and returns a
  [bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  object, which has [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and other methods.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  returns sampler and grid diagnostics.

- [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md)
  summarizes the estimated distribution of effects, and
  [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  and
  [`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  summarize the site effects.

- [plot()](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  plots a fit. The functions that it calls can also be called directly:
  [`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
  [`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md),
  [`plot_shrinkage()`](https://joonho112.github.io/bayesEfron/reference/plot_shrinkage.md),
  [`plot_exceedance()`](https://joonho112.github.io/bayesEfron/reference/plot_exceedance.md),
  [`plot_site_posterior()`](https://joonho112.github.io/bayesEfron/reference/plot_site_posterior.md),
  [`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md),
  [`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md),
  [`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md),
  [`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md)
  and, for simulated data,
  [`plot_recovery()`](https://joonho112.github.io/bayesEfron/reference/plot_recovery.md).

- [predict()](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md)
  returns draws of the effect or of the estimate of a new site.

- [loo()](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  estimates how well a fit predicts the estimate of a site that is left
  out, and
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  compares fits on several grids by that measure.

- [`as_normalmix()`](https://joonho112.github.io/bayesEfron/reference/as_normalmix.md)
  returns the estimate of \\g\\ as a `normalmix` object for the ebnm
  package.

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  compiles the two Stan programs before the first fit, and
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  deletes the compiled models.

## Data sets and vignettes

The following data sets come with a stored fit, which can be used
without CmdStan.

- [krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md)
  and
  [krw_fit](https://joonho112.github.io/bayesEfron/reference/krw_fit.md):
  the callback gaps of 97 firms, analyzed in
  [`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md).

- [raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
  and
  [raudenbush_fit](https://joonho112.github.io/bayesEfron/reference/raudenbush_fit.md):
  a meta-analysis of 19 studies, analyzed in
  [`vignette("teacher-expectancy")`](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md).

- [twin_towers](https://joonho112.github.io/bayesEfron/reference/twin_towers.md)
  and
  [twin_towers_fit](https://joonho112.github.io/bayesEfron/reference/twin_towers_fit.md):
  a simulated data set with known true effects, analyzed in
  [`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md).

The data set
[star_schools](https://joonho112.github.io/bayesEfron/reference/star_schools.md)
contains estimates of the effect of small classes in kindergarten in 79
schools of Project STAR, a multisite trial.

The other vignettes cover data input, the choice of the grid,
diagnostics, plots, a comparison with empirical Bayes methods, the
model, and the implementation. `browseVignettes("bayesEfron")` lists
them all.

## Funding

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama. The opinions expressed are those of the authors and do not
represent views of the Institute or the U.S. Department of Education.

## See also

Useful links:

- <https://github.com/joonho112/bayesEfron>

- <https://joonho112.github.io/bayesEfron/>

- Report bugs at <https://github.com/joonho112/bayesEfron/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
