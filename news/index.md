# Changelog

## bayesEfron 0.4.0

The defaults of
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
have changed, so the same data and seed no longer give the posterior
draws of 0.3.0. By default, a short preliminary fit (the pilot fit) now
chooses the range of the grid, a new Stan program samples the model, and
a fit stores the draws of the distribution of effects but no draws for
each site. Use `settings = "0.3.0"` to fit the model as 0.3.0 did, with
the same draws if the version of CmdStan and the machine are the same. A
fit saved by an earlier version can still be used: its summaries are
computed again from its draws with the definitions of this version, and
a warning names what was recomputed.

### Breaking changes

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  and
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  have a new default grid method, `grid_method = "pilot"`: a pilot fit
  on the range of 0.3.0 narrows that range to where the distribution of
  effects lies (see
  [`?make_efron_grid`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)).
  With this default,
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  needs CmdStan. Use `grid_method = "paper_realdata"` for the range of
  0.3.0 without a fit. The first fit of a session that relies on the
  default says in a message that the defaults have changed.
- With the default `L = NULL`,
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  and
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  now choose the number of grid points from the data: 101, or more (up
  to
  401. when 101 are too few for the smallest standard error or for the
       spread of the effects. The default was 101 points for every data
       set. The largest `L` that can be given is now 401 (was 300).
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  runs a new Stan program by default. For the same seed its draws differ
  from the draws of 0.3.0.
- A fit no longer stores five variables per site (`theta_mean`,
  `theta_sd`, `theta_map`, `theta_rep` and `log_lik`);
  `site_draws = TRUE` stores them. The methods of a fit, `loo()` and
  `waic()` included, compute what they need from the draws of the
  distribution of effects.
- [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) now
  compute the limits of a site’s interval as quantiles of the posterior
  distribution of its effect, which averages over the draws of the
  distribution of effects. The limits were sample quantiles of the
  `theta_rep` draws. Those draws take values on grid points, although
  R’s interpolation can put their sample quantiles between points.
  `confint(method = "draws")` gives those limits for a fit made with
  `site_draws = TRUE`.
- With point masses on the grid, quantiles and tail probabilities are
  now computed with the probability of each grid point spread evenly
  over one grid spacing around it (`reading = "interpolated"`).
  [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md),
  [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
  [`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  and [`quantile()`](https://rdrr.io/r/stats/quantile.html) have a
  `reading` argument, and `reading = "step"` leaves the probability on
  the grid point.
- In the table of site effects, the columns `hpdi_lower`, `hpdi_upper`
  and `map` are now `lower`, `upper` and `mode`, and the table gains the
  columns `label`, `estimate` and `se`.
  [`confint()`](https://rdrr.io/r/stats/confint.html) returns `mean` in
  place of `point` and, for `type = "g"`, `parameter` in place of
  `site`. The earlier columns remain with the same values and will be
  removed in 0.5.
- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and the tables now
  name the sites by their labels when the data had labels (the sites
  were numbered before). In the `parm` argument of
  [`confint()`](https://rdrr.io/r/stats/confint.html), numbers are
  positions and character strings are labels; `parm = "3"` was read as
  the third site.
  [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  and
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  now throw an error for duplicated site labels.
- The plots use ggplot2, which the package now imports. They carry a
  caption in place of a title, name the sites by their labels and list
  them from the top. `plot(type = "g")` now shows the density of the
  distribution of effects in place of the probabilities of the grid
  points, which `plot_g(type = "mass")` shows.
- `plot(backend = "base")` is deprecated and will be removed in 0.5.
  Until then it plots the four types of 0.3.0 with base graphics and
  issues a warning once per session.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) also no
  longer attaches the plotted numbers to its result (the undocumented
  attribute `bef_plot_payload`);
  [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md),
  [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
  [`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  and
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  return them.
- [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) now throw an
  error for an argument that they do not have (it was ignored).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) still accepts
  `level` and `sort_by` for the four types of 0.3.0.
- [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) now throw an error that
  points to `loo()`. The numbers that they returned were computed from a
  posterior mean of the log likelihood and are not information criteria.
- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  now warns about an E-BFMI below 0.3, the value below which CmdStan
  warns (the threshold was 0.2). An infinite R-hat now fails the check
  instead of putting R-hat among the diagnostics that could not be
  computed. A quantity has an infinite R-hat when it stays at one value
  in each chain and the values differ between chains.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  now shows the output of CmdStan and the cmdstanr messages about the
  sampler in an interactive session only. Use `show_messages = TRUE` to
  see them elsewhere.
- In
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md),
  the arguments after `theta_hat` and `sigma` must now be named. The
  grid method `"kl_target_experimental"` and the `kappa` argument are
  removed, and a grid no longer has a `kappa` element.
- A result of [`summary()`](https://rdrr.io/r/base/summary.html) that
  was saved by an earlier version can no longer be printed. Call
  [`summary()`](https://rdrr.io/r/base/summary.html) again on the fit.

### New features

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  has a new `settings` argument: `settings = "0.3.0"` runs the Stan
  program of 0.3.0 with the range, the 101 grid points and the
  site-level draws of that version. An argument given with it keeps its
  value.
  [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  compiles both Stan programs, and
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  deletes both compiled models.
- Three new functions summarize a fit.
  [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md)
  returns the density, the distribution function, quantiles and tail
  probabilities of the distribution of effects, with intervals.
  [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  returns the table of site effects at any level, and
  [`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
  the probability that a site effect exceeds a threshold.
- Ten new functions plot a fit:
  [`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md),
  [`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md),
  [`plot_shrinkage()`](https://joonho112.github.io/bayesEfron/reference/plot_shrinkage.md),
  [`plot_exceedance()`](https://joonho112.github.io/bayesEfron/reference/plot_exceedance.md),
  [`plot_site_posterior()`](https://joonho112.github.io/bayesEfron/reference/plot_site_posterior.md),
  [`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md),
  [`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md),
  [`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md),
  [`plot_recovery()`](https://joonho112.github.io/bayesEfron/reference/plot_recovery.md)
  and
  [`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md).
  They return ggplot objects, except that `plot_mcmc(type = "pairs")`
  returns a grid of plots.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a fit
  calls them through `type`, where `"sites"` and `"diagnostics"` replace
  `"caterpillar"` and `"diagnostic"`; the earlier names are still
  accepted.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  and
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  accept a data frame or a `bef_data` object as their first argument.
- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  takes the labels of a metafor `escalc` object from its `slab`, and its
  result has an
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method.
  For a data frame with variances and no standard errors it takes the
  square roots, with a message. It warns when a column with a name that
  begins with `vi` or `var` is given as the standard errors.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  and
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  have a new `margin` argument for the `"pilot"` method, which sets how
  far the region found by the pilot fit is extended on each side, as a
  fraction of its width.
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  also has a `seed` argument for the pilot fit.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  and
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  gain a `support` argument, which gives the two limits of the grid
  directly, and a `knots` argument, for which `"grid"` is the only value
  so far.
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  also takes a grid made by
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  as `grid`. A grid has the class `bef_grid`, with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
- `component = "normal"` puts a normal distribution with a standard
  deviation of `width` grid spacings at each grid point, in place of a
  point mass, as in the penalized Gaussian mixture of Sui and Tipton
  (2026).
- A warning of class `bef_warn_grid_resolution` reports a grid that is
  too coarse for the smallest standard error or for the spread of the
  effects, or narrower than the median standard error.
- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  reports, for each end of the grid, the largest posterior probability
  that a site has in the outer 2% of the range, and the posterior
  distribution of the standard deviation of the effects together with
  its prior.
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  issues a warning of class `bef_warn_grid_ends` when that probability
  exceeds 0.05 at an end that was specified in `support` or that follows
  from the estimates alone. The result of
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  has an [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  method.
- [`confint()`](https://rdrr.io/r/stats/confint.html) has the new
  arguments `method` and `interval`, with `interval = "hpd"` for the
  shortest interval, and its result prints the level and the kind of
  interval.
- A fit has the new methods
  [`quantile()`](https://rdrr.io/r/stats/quantile.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) for joint draws
  of the site effects,
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`posterior::ndraws()`](https://mc-stan.org/posterior/reference/draws-index.html)
  and
  [`posterior::variables()`](https://mc-stan.org/posterior/reference/variables.html).
  `vcov(full = TRUE)` returns the covariances between the site effects,
  and `coef(type = "mode")` is the new name of `type = "map"`, which is
  still accepted.
- [`print()`](https://rdrr.io/r/base/print.html) of a fit shows the
  number of draws, the range of the standard errors, the grid with the
  method that chose its range, and the functions to call next.
- [`predict()`](https://rdrr.io/r/stats/predict.html) accepts a data
  frame of estimates and standard errors as `newdata`. It then returns
  the posterior summaries of new sites with known estimates, with
  intervals at `level`: for each row, the summaries that the site would
  have in a fit of the data with that site added.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  gains the sampler arguments `refresh`, `init`, `thin` and
  `show_messages`. `thin` reduces what is stored; the sampler
  diagnostics still describe the whole run.
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  accepts several grid methods, ranges given in `support` and both kinds
  of component, and `keep_fits = TRUE` returns the fits as well. For
  each specification it also returns the range, the mean and the
  standard deviation of the distribution of effects, the posterior means
  of the sites named in `sites`, the posterior probability at the ends
  of the grid (as in
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)),
  and whether the grid is fine enough. It reports which fit is running,
  and its result has a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method. A
  comparison saved by 0.3.0 cannot be plotted and has to be run again.
- New
  [`as_normalmix()`](https://joonho112.github.io/bayesEfron/reference/as_normalmix.md)
  returns the posterior mean of the distribution of effects as a
  `normalmix` object, which the ebnm package accepts as a fixed prior.
  With the bayesplot package installed, `nuts_params()`,
  `log_posterior()`, `rhat()` and `neff_ratio()` work on a fit. Both
  packages are suggested: bayesplot for the plots of the sampler, and
  ebnm for the comparison with empirical Bayes methods.
- `krw_discrimination` and `twin_towers` are new data sets.
  `krw_discrimination` contains the firm-level callback gaps of Kline,
  Rose and Walters, from the public replication files of Walters (2024),
  as analyzed by Lee and Sui (2025). Both data sets record where they
  come from in the attribute `source`, and the scripts that rebuild them
  from the public files are in `data-raw/` in the source repository.
- `star_schools` is a new data set from a multisite trial in education:
  the estimated effect of small classes in each of the 79 kindergarten
  schools of Tennessee’s Project STAR. The estimates are built from the
  public student file with a model that has a random intercept for each
  class. The data set also has the estimates and standard errors that
  take the students of a school as independent.
- The stored fits `krw_fit` and `twin_towers_fit` are new. Like
  `raudenbush_fit`, they are made with the grid method, the Stan program
  and the storage of this version (`twin_towers_fit` with 8 spline
  functions) and are thinned to 800 of 12,000 draws (`thin = 15`), so
  that the examples run without CmdStan. `raudenbush_fit` has the range
  that it had in 0.3.0 but different draws.

### Bug fixes

- `predict(type = "theta_hat")` with several standard errors now gives
  each new site an effect of its own, drawn from the same posterior draw
  of the distribution of effects. All the new sites of a draw had the
  same true effect, so that their estimates differed by sampling error
  alone.
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  no longer fails with an internal error when one of its fits fails. It
  issues a warning that names the specification and returns a row of
  missing values for it, as documented.
- The errors of
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  and
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  about the data now state in one sentence what is wrong: a missing
  value, an infinite value, a length or a sign. A missing value was
  reported inside a sentence about the type and the length of the
  vector.
- `make_efron_grid(grid = )` no longer throws an error about
  `grid_method`, the argument that the name was matched to. Its error
  now points to `support` and `L`, and `compare_efron_grids(grid = )`
  reports that a grid cannot be given.
- In the posterior standard deviation of a site, the variance between
  draws is now divided by the number of draws, as its definition
  requires (was the number of draws minus one).

### Documentation

- The help pages and the vignettes are rewritten in plainer language.
  All ten vignettes are built around the data sets and stored fits
  listed above, with the callback gaps (`krw_discrimination`) as the
  main example. Four are new: `twin-towers`, on the recovery of a
  distribution with two groups; `empirical-bayes`, a comparison with the
  ebnm and deconvolveR packages; `plots`, a gallery of the plot
  functions; and `teacher-expectancy`, which keeps the 19 studies that
  were the main example of 0.3.0 as a second case, with a simulation of
  what the model can and cannot do with so few.
- The article “Small classes in 79 schools” on the package website
  analyzes `star_schools`: how the effect varies across schools
  according to the fit, how much of that rests on the standard errors,
  and what simulated data sets with known effects show. The article and
  the fit that it uses are not installed with the package.
- The installed files `benchmark-coverage.csv` and
  `benchmark-remedies.csv` (in `verification/`) and
  `grid_comparison_example.rds` (in `examples/`) are replaced by the
  tables and the stored comparisons that the vignettes of this version
  use.

## bayesEfron 0.3.0

This is the first public release since 0.1.0. It includes the changes of
0.2.0, which was not released separately and is described below. The
model is the same as in 0.2.0: with the same data, settings, seed and
version of CmdStan on the same machine, the posterior draws are
identical.

### Breaking changes

- In the table of site effects that
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return,
  `sd` is now the posterior standard deviation of each effect and `map`
  its posterior mode, as defined in Lee and Sui (2025, Section 4.3 and
  Appendix A). Earlier versions gave the standard deviation and the mode
  for a given distribution of effects, averaged over the posterior draws
  of that distribution. That `sd` left out the uncertainty about the
  distribution (the posterior standard deviations of the example data
  are 3% to 35% larger), and that `map` was not, in general, the mode of
  the posterior distribution.
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and `coef(type = "map")`
  change with them. The posterior means and the intervals are unchanged.
- A fit saved by an earlier version holds `sd` and `map` in the earlier
  definitions, and one saved by 0.1.0 holds `effective_params` in its
  first definition too.
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and the other methods of
  a fit compute these again from the stored draws and give a warning
  that names what changed; the saved object is not altered. A saved
  result of [`summary()`](https://rdrr.io/r/base/summary.html) cannot be
  brought up to date and has to be made again from the fit. A fit saved
  by 0.1.0 still has to be fitted again before
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
  `loo()`, `waic()` or
  [`predict()`](https://rdrr.io/r/stats/predict.html) can be used.
- The compiled model is now handled by cmdstanr alone. It is stored in
  `tools::R_user_dir("bayesEfron", "cache")` under a name made of a hash
  of the Stan file, the CmdStan version and the machine architecture,
  and is recompiled when one of them changes. The lock files and the
  files with metadata that earlier versions wrote next to it are gone.
  The environment variable `BAYESEFRON_CACHE_ROOT` still moves the
  directory.
- [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  deletes the compiled models and returns their number. The values
  `"session"`, `"compiled_models"` and `"lock_only"` of its argument
  `scope` are accepted with a warning and will be removed; `"lock_only"`
  no longer does anything. Files left by earlier versions are deleted as
  well.
- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  no longer has the argument `seed_for_check`.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) no longer has
  `type = "sensitivity"`, which drew the same picture as
  `type = "caterpillar"`.
- Three environment variables are replaced. Use `plot(backend = "base")`
  for `BAYESEFRON_NO_GGPLOT2`, `options(bayesEfron.use_cli = FALSE)` for
  `BAYESEFRON_NO_CLI`, and the argument `parallel_chains` of
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  for `BAYESEFRON_PARALLEL_CHAINS`.
- The environment variables that switched parts of the test suite on
  (`BAYESEFRON_RUN_LIVE`, `BAYESEFRON_RUN_FULL_LIVE`,
  `BAYESEFRON_RUN_PARITY`, `BAYESEFRON_TIER3_OK_TO_REFIT`,
  `BAYESEFRON_TIER3_FULL_LIVE_MATRIX`, `BAYESEFRON_SOURCE_ROOT` and
  `BAYESEFRON_ENFORCE_RELEASE_LEDGER`) are no longer read. The tests
  that fit a model run when CmdStan is installed and are skipped when it
  is not.
- Condition classes that were never signaled are removed, as are the
  elements `stage`, `module` and `predicate` of error objects. The class
  `bayesEfron_validate_error` is now `bef_validate_error`. Every error
  of the package still inherits from `bef_error`.
- The element `attribution` of the list that
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  returns names the section of Lee and Sui (2025) that a rule comes
  from.
- The packages ps, jsonlite, deconvolveR and desc are no longer
  suggested.

### New features

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  has the arguments `max_treedepth` and `parallel_chains`.
- A fit records the settings of the sampler and the versions of
  bayesEfron and cmdstanr that made it, as attributes of `metadata` (see
  [`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) has the
  argument `backend`, `"ggplot2"` or `"base"`.
- The data set `raudenbush1985`, 19 studies of teacher expectancy, and
  `raudenbush_fit`, the model fitted to it with the default settings.
  All examples and vignettes use them and run without CmdStan.
- [`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  describes the fitted object and lists its methods.
- A compiled model that does not run is detected and compiled again (on
  Windows only an empty file is detected).

### Bug fixes

- `predict(type = "theta_hat", n = 1)` stopped with an error when
  several standard errors were given. It returns a matrix with one row.
- The variance of an effect given the distribution of effects, from
  which `sd` and `effective_params` are computed, was taken as a
  difference of two moments. That loses its digits when the estimates
  are very far from zero compared with their standard errors. The
  variance is now taken about the mean. Results for estimates of
  ordinary size are unchanged.
- `loo()` returned `NaN` as the Monte Carlo standard error of a site
  whose likelihood was so small that it rounded to zero.
- A fit with fewer than 400 draws was not checked at all. It is now
  checked for divergent transitions, iterations at the maximum tree
  depth and E-BFMI. R-hat and the effective sample sizes are still not
  compared with their thresholds in such a fit, and `diagnostic_skipped`
  says so.
- A fit in which R-hat or an effective sample size was infinite for some
  quantity stopped with an error. The diagnostic is now stored as `NA`
  and named in `diagnostic_skipped`.
- `iter_sampling = 0` was accepted, and the function stopped with an
  error only after the warmup had run. It is now rejected with the other
  arguments.
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  put `elpd_diff` and `se_diff` in the wrong rows with loo 2.10 or
  later, so that the best specification could show the largest
  difference. They are now computed from the pointwise values.
- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  called [`set.seed()`](https://rdrr.io/r/base/Random.html), which reset
  R’s random number generator for the rest of the session and made two
  fits with the same seed write to the same temporary files. It no
  longer does.
- `plot(type = "g")` drew the summaries of the distribution, not the
  distribution, for a fit made with the default settings. It now
  computes the probabilities from the draws.
- `plot(type = "diagnostic")` with ggplot2 drew five diagnostics of
  different scales on one axis. Each has its own panel.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) showed the names of
  the same diagnostics under three labels. There are now two, failed
  checks and warnings.
- Errors named an internal helper as their source. They name the
  function that was called.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) stopped when
  given further arguments; it ignores them, as the other methods do.
- `predict(seed = )` left R’s random number generator seeded when the
  session had not used it before. It now leaves the generator as it
  found it.
- `plot(type = "loo", backend = "base")` drew the reference lines at 0.5
  and 0.7 only when the values reached them. Both backends now show them
  always.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) showed the number
  of iterations at the maximum tree depth under the label
  `max treedepth`, which read as a depth. The label is now
  `treedepth hits`. The titles and axis labels of the plots say what is
  drawn (“Distribution of effects” and so on).
- The stored test inputs held the local paths of the files they were
  made from.

### Documentation

- The help pages are rewritten. They say what the columns of the table
  of site effects are. In particular `hpdi_lower` and `hpdi_upper` are
  the 5% and 95% quantiles of the posterior draws and not the limits of
  a highest posterior density interval.
- The README says how to install the vignettes with the package; an
  installation from GitHub does not build them unless asked to.
- Six vignettes replace sixteen: `bayesEfron`, `data-input`,
  `choosing-a-grid`, `diagnostics`, `model` and `implementation`.
- [`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)
  shows how far the results for the example data depend on the grid.
  With 19 studies the posterior mean of the study with the largest
  estimate changes with the width of the grid and the number of spline
  functions, while the summaries of the distribution of effects and the
  other studies hardly change.
- The table of coverage on the benchmark of Lee and Sui (2025) reports,
  for each number of sites, how many fits did not converge. With 1,500
  sites and the default settings three of twenty did not, and a second
  table shows what a narrower grid and a larger maximum tree depth
  change for those three. The scripts and the data of the benchmark are
  in `data-raw/benchmark/` in the source repository.

## bayesEfron 0.2.0

Not released separately.

### Breaking changes

- `effective_params` is the sum over sites of the posterior variance of
  the effect, for a given distribution of effects, divided by the
  squared standard error. Version 0.1.0 reported the number of sites
  minus that sum.
- `fit$posterior` holds only the draws of `mean_g`, `var_g`, `sd_g`,
  `effective_params` and `log_marginal_likelihood`. The draws for each
  site are in `fit$draws`, and `fit$metadata$theta_rep_draws` is
  removed.
- The draws no longer include `g` and `log_w`;
  `store_grid_quantities = TRUE` puts them back.
- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  returns R-hat and the effective sample sizes for every quantity in the
  draws, where 0.1.0 returned one number for all of them together, and
  it returns E-BFMI. R-hat above 1.05 or a divergent transition is a
  failed check; R-hat above 1.01, an effective sample size below 400, an
  iteration at the maximum tree depth and E-BFMI below 0.2 are warnings.
- Objects saved by 0.1.0 can still be used with
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
  [`print()`](https://rdrr.io/r/base/print.html).
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
  `loo()`, `waic()` and
  [`predict()`](https://rdrr.io/r/stats/predict.html) need the model to
  be fitted again.

### New features

- `loo()` and `waic()` methods, and `plot(type = "loo")`.
- [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  fits several specifications of the grid and of the spline basis and
  ranks them by leave-one-out cross-validation.
- [`predict()`](https://rdrr.io/r/stats/predict.html) draws the effect
  or the estimate of a new site, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
  posterior means.
- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  accepts a data frame, finds the columns by name and takes variances as
  well as standard errors. It stops on missing values.

### Performance

- The sampler is about 12 to 34 times faster than in 0.1.0 for 200 to
  1,500 sites, and a fitted object is about half the size. The normal
  densities of the estimates at the grid points are computed once, and
  the likelihood is one matrix product.

## bayesEfron 0.1.0

First public release.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  fits the fully Bayesian Efron log-spline model to estimates and their
  standard errors.
- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  builds the grid and the spline basis, with three rules from Lee and
  Sui (2025) and one experimental rule.
- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  converts a list or a
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  object.
- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  and
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  manage the compiled model.
- Methods for fitted objects:
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html),
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
