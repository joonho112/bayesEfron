# Changelog

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
