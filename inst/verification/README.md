# Simulated data sets with known effects

Three tables, computed with version 0.4.0. `vignette("implementation")`
sets them side by side.

## Independent panels (version 0.4.0)

Every panel has its own true effects, drawn from the two uniform
populations of the Twin Towers example (a third of the sites uniform on
-1.7 to -0.7, two thirds uniform on 0.7 to 2.7), with reliability 0.7 and
the largest standard error 9 times the smallest. There are 20 panels of
1,500 sites and 100 panels of 100 sites. All intervals are equal-tailed
90% intervals.

`panels-convergence.csv` has one row for each of three settings of
`bayes_efron_fit()` on the 20 panels of 1,500 sites (`setting`): the
defaults, 8 spline functions (`M = 8`), and `settings = "0.3.0"`.
`panels` is the number of panels. `valid` counts the fits that meet three
conditions: the largest R-hat over `alpha`, `lambda`, `mean_g` and `sd_g`
is below 1.01, the smallest bulk effective sample size is at least 400,
and there is no divergent transition. `rhat_max`, `rhat_above_1.01` and
`rhat_above_1.05` describe that R-hat over the panels, `ess_bulk_min` is
the smallest bulk effective sample size, and `divergences` is a total.
`treedepth_share_median` and `treedepth_share_max` are the proportions of
a fit's 12,000 iterations that reached the maximum tree depth,
`leapfrog_median` the number of leapfrog steps of a fit, and
`sampling_seconds_median` and `pilot_seconds_median` the times of
sampling and of the short pilot fit of the `"pilot"` method. `rmse`,
`inclusion` (the proportion of the true effects inside the intervals),
`mean_length` and `w1_population` (the Wasserstein distance between the
estimate of g and the population of the true effects) are means over the
panels, so `inclusion` estimates the coverage of the intervals.
`inclusion_se` is the standard deviation of that proportion over the
panels divided by the square root of their number.

`panels-methods.csv` sets the fit with the defaults beside three other
ways of making intervals from the same data (`method`), for both numbers
of sites (`K`, with `panels` the number of panels):
the posterior mean of g taken as known (`g of bayesEfron, plugged in`),
and the empirical Bayes fits `ebnm::ebnm_normal(mode = "estimate")` and
`ebnm::ebnm_npmle()`. Each of `rmse`, `inclusion` and `mean_length` is a
mean over the panels, followed by its Monte Carlo standard error.

In all three tables of this directory a column with a name ending in
`_se` is the Monte Carlo standard error of the column before it,
`generated_on` is the date on which the table was computed and
`cmdstan_version` the version of CmdStan that was used.

The script is `data-raw/ebnm-comparison.R` in the source repository; the
fits of the panels are made by `data-raw/benchmark/panels-run.R`. The
tables are not recomputed when the package is checked, because the fits
take about two hours.

## Nineteen studies with known effects (version 0.4.0)

`few-studies.csv` is the table of `vignette("teacher-expectancy")`. Every
data set has the 19 standard errors of `raudenbush1985`, and an estimate
is its true effect plus normal noise with its standard error. The true
effects are drawn from one of three distributions (`scenario`): a normal
distribution with the mean and the standard deviation that REML estimates
for the 19 studies, a normal distribution with the same mean and the
posterior mean of the standard deviation of g in `raudenbush_fit` as its
standard deviation, and the posterior mean of g in that fit. `true_sd` is
the standard deviation of the distribution. There are 200 data sets for
each (`replications`).

Each data set is analyzed in four ways (`method`): `bayes_efron_fit()`
with its defaults; the same with `adapt_delta = 0.99`; the normal model
with its mean and standard deviation fixed at their REML estimates; and
the posterior distribution under the distribution that the effects were
drawn from (`known g`).

`sd_mean`, `sd_se`, `sd_q05`, `sd_q50` and `sd_q95` describe the estimate
of the standard deviation of the effects over the data sets (the
posterior mean of the standard deviation of g, or the REML estimate),
`sd_zero` is the proportion of data sets in which it is zero, and
`sd_inclusion` the proportion in which the 90% interval for the standard
deviation of g contains the true value. `inclusion`, `length` and `rmse`
are the proportion of the true effects inside the 90% intervals, the mean
length of those intervals and the root mean squared error of the
posterior means, each a mean over the data sets followed by its Monte
Carlo standard error. `largest_inclusion`, `largest_error` and
`largest_error_se` concern the study with the largest estimate in each
data set: the proportion of data sets in which its interval contains its
true effect, and the mean of its posterior mean minus its true effect.

For the fits of this package, `rhat_above_1.01`, `with_divergences` and
`with_many_divergences` count the fits with an R-hat above 1.01, with at
least one divergent transition and with more than 10;
`divergences_median` and `divergences_max` describe the number of
divergent transitions of a fit; `inclusion_with_divergences` and
`inclusion_without_divergences` give the coverage in the fits with and
without them; `end_warnings` counts the fits with a warning about an end
of the grid; and `width_mean` is the mean width of the grid. In the rows
for `adapt_delta = 0.99`, `sd_change`, `inclusion_change` and
`length_change` are mean differences from the fit with the defaults on
the same data set, each followed by its standard error.

`observed_own_sd` and `observed_reml_sd` are the two estimates of the
standard deviation for the 19 studies themselves. `own_as_large` is the
proportion of data sets of the scenario in which the estimate of this
package is at least as large as the first, `reml_as_small` the proportion
in which the REML estimate is at most the second, and `both_as_observed`
the proportion in which both conditions are met.

The script is `data-raw/benchmark/few-studies.R` in the source
repository. It makes 1,200 fits and takes about ten minutes on 14 cores;
the table is not recomputed when the package is checked.

For the third distribution an effect is a grid point drawn with its
probability and moved uniformly within its cell of the grid. Its
`known g` row therefore uses a prior distribution that is uniform within
each cell: the normal likelihood is averaged over each cell, and within a
cell the posterior distribution is a normal distribution truncated at the
limits of the cell.

The results of every data set, four analyses for each of the 600, are in
`data-raw/benchmark/results-v040/few-studies-rows.csv` in the source
repository. The `known g` rows need no fit. For the third distribution,
`Rscript data-raw/benchmark/few-studies.R --analytic-only` computes them
and their row of this table again, and
`data-raw/benchmark/histogram-posterior-check.R` compares the formulas
with numerical integration.
