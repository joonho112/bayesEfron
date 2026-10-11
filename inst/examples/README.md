# Stored results that the documentation reads

This directory contains four files. None of them is needed to use the
package. The vignettes and the examples in the help pages read them,
together with the three fits in `data/`, so that no model is fitted and no
other package is needed when the documentation is built.
`twin_towers_comparison.rds` is made by `data-raw/example-fits.R` together
with the fits, `lee_sui_examples.rds` by `data-raw/example-evidence.R`,
`raudenbush_ranges.rds` by `data-raw/range-comparison.R`, and
`ebnm_comparison.rds` by `data-raw/ebnm-comparison.R`; the scripts are in
the source repository.

## The stored fits

| Object | Sites | Call | Range of the grid | Draws run / stored |
|:--|--:|:--|:--|:--|
| `krw_fit` | 97 | `bayes_efron_fit(as_bef_data(krw_discrimination, site = "firm_name"), seed = 20250807, thin = 15)` | -0.083 to 0.159, the range of `"paper_realdata"` | 12,000 / 800 |
| `raudenbush_fit` | 19 | `bayes_efron_fit(raudenbush1985, seed = 1985, thin = 15)` | -1.07 to 1.93, the range of `"paper_realdata"` | 12,000 / 800 |
| `twin_towers_fit` | 1,500 | `bayes_efron_fit(twin_towers, M = 8, seed = 20250627, thin = 15)` | -3.00 to 3.67, 30% of the range of `"paper_realdata"` | 12,000 / 800 |

All three were fitted with CmdStan 2.38.0 on 2026-10-10 with the defaults
of this version: 4 chains of 1,000 warmup and 3,000 sampling iterations,
the `"pilot"` method for the range, 101 grid points. The fits store the
draws of g and no site-level draws. The diagnostics that they store
describe all 12,000 draws: the largest R-hat is 1.003, the smallest bulk
effective sample size is 994, and there is no divergent transition and no
iteration at the maximum tree depth. The summaries are computed from the
800 stored draws. Compared with the summaries from all 12,000 draws, a
limit of a 90% site interval moves by at most 0.44, 0.58 and 0.27 grid
spacings in the three fits, and a site's probability of a positive effect
by at most 0.004.

## `twin_towers_comparison.rds`

A list. Its element `table` has one row for each of five fits of the 1,500
simulated sites at full length (12,000 draws, seed 20250627): the call of
`twin_towers_fit`; the defaults, with 6 spline functions; the range
around the true effects (`"paper_simulation"`), with 6 and with 8 spline
functions; and the range of `"paper_realdata"`, which uses no true
effects and is not narrowed. The columns give the range, the
diagnostics, the time that sampling took, and how well the fit recovers
the true effects: `w1`, the Wasserstein distance between the estimate of g
and the true effects; `rmse`, the root mean squared error of the posterior
means; and `coverage_90`, the proportion of the true effects inside the
90% intervals. `w1_over_reference` and `rmse_over_reference` divide by the
values of the fit named in `reference`, the range around the true effects
with 6 spline functions. This is one simulated panel; the numbers do not
establish coverage over simulated data sets.

## `lee_sui_examples.rds`

A list with five elements.

- `sensitivity`: one row for each of six fits at full length. For the
  callback gaps: the default fit (`baseline`), a range that extends beyond
  the estimates by a quarter of their width (`narrower`), and 8 spline
  functions (`basis8`). For the simulated sites: the rows `baseline`,
  `default_M6` and `true_range` of the comparison above.
- `settings`: the settings and hashes of four of those fits: the three
  fits of the callback gaps and the `baseline` fit of the simulated sites.
- `curves`: the estimate of g on its grid for the three fits of the
  callback gaps and for `twin_towers_fit`, as probabilities and as
  probabilities divided by the grid spacing.
- `cases`: for each of the two data sets, the normal random-effects model
  fitted by REML with `metafor::rma()` to the same estimates, with plug-in
  site intervals, which fix its mean and variance at their estimates; the
  probability of a positive effect for every site from the stored fit and
  from the full run, both with `reading = "step"`, and from the stored fit
  with the default reading; and summaries of both.
- `about`: dates, versions and the names of the scripts.

## `raudenbush_ranges.rds`

The result of `compare_efron_grids()` for the 19 teacher expectancy
studies on three ranges: the range of the estimates widened on each side
by 0.25, 0.5 and 1 times its width, with 101 grid points and 6 spline
functions (seed 1985, 12,000 draws for each fit). It stores the estimate
of g from every fit and the summaries of studies 4 and 8, so that `plot()`
can show the comparison. The examples of `?compare_efron_grids` and
`?plot.bef_grid_comparison` read it.

## `ebnm_comparison.rds`

A comparison of the stored fits with empirical Bayes fits of the same
data. All intervals are equal-tailed 90% intervals. For the ebnm fits they
are the sample quantiles of 4,000 draws from each posterior distribution.
A list:

- `krw`: for the 97 callback gaps, `table` has one row for each of six
  empirical Bayes methods of the ebnm package, for the posterior mean of g
  of `krw_fit` taken as known, and for `krw_fit`. Its columns give the
  correlation of the posterior means with the posterior means of
  `krw_fit` and their largest difference in units of the site's standard
  error, the mean posterior standard deviation and the mean length of the
  intervals with their ratios to the values for `krw_fit`, and the
  maximized log likelihood. `firms` has the 10 firms with the largest
  posterior means under each method, and `cdf` the distribution function
  of each estimate of g, with the band of `krw_fit`.
- `twin`: for the 1,500 simulated sites, `table` adds the true
  distribution taken as known and a fit on the range around the true
  effects, and gives the root mean squared error of the posterior means
  (over all sites and within the thirds of the standard errors), the
  proportion of the true effects inside the intervals with its standard
  error, the mean length, the Wasserstein distance between the estimate of
  g and the true effects, and the log likelihood. `cdf` is as above.
- `equal_se`: the same true effects with standard errors of 1. `table`
  compares the fit with `deconvolveR::deconv()`, run on the fit's grid and
  with its number of spline functions; `g` has the two estimates of g with
  their bands; and `agreement` shows how close they are. `repeated` gives
  the same measures over 20 such data sets, with 6 and with 8 spline
  functions.
- `panels`: `table` gives, for 100 panels of 100 sites and 20 panels of
  1,500 sites, the error, the coverage of the intervals and their mean
  length for the fit, its g taken as known, and the normal and
  nonparametric empirical Bayes fits, each with its Monte Carlo standard
  error; `by_panel` has the values of every panel; and `convergence` the
  diagnostics of the fits of the large panels under three settings.
- `settings`: the call of every method.
- `about`: the seeds, the versions of the packages and the names of the
  protocol and the script.

A flat prior has no distribution function and no log likelihood, so those
entries of its rows are missing. These are single data sets and a limited
number of panels from one data-generating process, so they do not
establish how the methods compare elsewhere.
