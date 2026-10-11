# Protocol of the comparison with empirical Bayes methods

Fixed on 2026-10-09, before the comparison was run. The methods, their
settings, the quantities and the measures below are not changed after
results are seen. Where one had to be changed, the change and the reason
are in the section "Changes" at the end.

The script is `data-raw/ebnm-comparison.R`. It needs the ebnm package
(1.1.42 here), deconvolveR and CmdStan. Its result is stored in
`inst/examples/ebnm_comparison.rds`, which the documentation reads, so
that the documentation is built without ebnm.

## Data

| Name | Sites | Use |
|:---|---:|:---|
| Callback gaps (`krw_discrimination`) | 97 | True effects unknown: agreement between the methods |
| Simulated sites (`twin_towers`) | 1,500 | True effects known: error, proportion of true effects inside the intervals, distance to the truth |
| Equal standard errors | 1,500 | The true effects of `twin_towers`, each with a new error of standard deviation 1 (seed 20261101): Efron's estimator under its own conditions |
| Independent panels (`data-raw/benchmark/panels_K1500.rds`) | 1,500 x 20 | Repeated: coverage with its Monte Carlo error |
| Independent panels (`data-raw/benchmark/panels_K100.rds`) | 100 x 100 | Repeated, few sites: where the uncertainty about g matters |

The panels come from `data-raw/benchmark/panels-make.R`: every panel has
its own true effects, uniform on (-1.7, -0.7) with probability 1/3 and
uniform on (0.7, 2.7) with probability 2/3; reliability 0.7, largest
standard error 9 times the smallest.

## Methods

| Label | Call | Notes |
|:---|:---|:---|
| No shrinkage | `ebnm_flat(x, s)` | |
| Normal | `ebnm_normal(x, s, mode = "estimate")` | The default mode of ebnm is 0 |
| Scale mixture of normals | `ebnm_normal_scale_mixture(x, s, mode = "estimate")` | |
| ash | `ebnm_ash(x, s, mode = "estimate", mixcompdist = "normal", prior = "uniform")` | The default of ashr penalizes the weights |
| Unimodal | `ebnm_unimodal(x, s, mode = "estimate")` | |
| NPMLE | `ebnm_npmle(x, s)` | Maximum over its own grid, which spans the range of the estimates |
| Known g (simulated data only) | `ebnm_npmle(x, s, g_init = <true g>, fix_g = TRUE)` | For `twin_towers`: point masses at its 1,500 true effects. For the panels: the population, as 2,000 equally weighted points at its quantiles |
| g of bayesEfron, plugged in | `ebnm_npmle(x, s, g_init = as_normalmix(fit), fix_g = TRUE)` | The posterior mean of g taken as known. Its posterior means, standard deviations and log likelihood are those of this call. Its intervals are computed from the same g by the package, which reads the point masses as it does for the fit (each spread over its cell), so that the two rows differ by the uncertainty about g and by nothing else |
| bayesEfron | `bayes_efron_fit()` | The stored fits `krw_fit` and `twin_towers_fit` (800 of 12,000 draws); for the panels, fits at full length |
| bayesEfron, range from the true effects (simulated sites only) | `bayes_efron_fit(twin_towers, grid_method = "paper_simulation", theta_true = , M = 8)` | A separate row, for reference; it uses what an analysis does not have |
| Efron's estimator (equal standard errors only) | `deconvolveR::deconv(tau, X, family = "Normal", deltaAt = NULL, pDegree = M, c0 = 1)` | On the grid of the bayesEfron fit of the same data and with the same number of spline functions |

All ebnm calls use `output = ebnm_output_all()`.

Not used: `ebnm_group()` (with unequal standard errors it uses the first
of each group), and `ebnm_deconvolver(x / s)` for data with unequal
standard errors (it estimates the distribution of the effects divided by
their standard errors, another target).

## Quantities

| Quantity | ebnm | bayesEfron |
|:---|:---|:---|
| Posterior mean | `coef()` | `site_summary()$mean` |
| Posterior standard deviation | `fitted()$sd` | `site_summary()$sd` |
| 90% interval, equal-tailed | `quantile(fit, probs = c(0.05, 0.95), nsim = 4000, seed = 1)`, the sample quantiles of 4,000 draws from the posterior (`confint.ebnm` is a 95% highest-density interval) | `site_summary()$lower`, `$upper` |
| Distribution function of g | From the fitted mixture, on a common grid of 401 points | `g_summary(type = "cdf")` with its band |
| Log likelihood | `logLik()`, the maximized marginal log likelihood | Through the plug-in row: the marginal log likelihood at the posterior mean of g. `logLik()` of a bayesEfron fit is a posterior mean, which is a different quantity |
| Time | Elapsed seconds of the call | Sampling seconds, and the seconds of the pilot fit |

## Measures

Callback gaps: correlation of the posterior means with those of
bayesEfron; largest absolute difference of the posterior means in units of
the standard error of the site; ratio of the mean posterior standard
deviation to that of bayesEfron; ratio of the mean length of the 90%
interval; plug-in log likelihood; time. For each method the 10 firms with
the largest posterior means.

Simulated sites: root mean squared error of the posterior means, over all
sites and within the thirds of the standard errors; proportion of the true
effects inside the 90% intervals, with its binomial standard error; mean
length of the intervals; Wasserstein distance between the estimate of g
and the empirical distribution of the true effects, on a common support;
plug-in log likelihood; time.

Equal standard errors: the same for Efron's estimator (plug-in posterior
of its estimate of g) and bayesEfron, and the largest difference between
the two estimates of g.

Panels: for each number of sites and each of bayesEfron, the g of
bayesEfron plugged in, Normal and NPMLE: the root mean squared error, the
proportion of the true effects inside the 90% intervals and the mean
length of those intervals, each as the mean over the panels with its Monte
Carlo standard error, the standard deviation over the panels divided by
the square root of their number.

## Large panels: convergence

The 20 independent panels of 1,500 sites were fitted with the
defaults, with `M = 8`, and with `settings = "0.3.0"`. The criterion: no
fit with the defaults whose largest R-hat over `alpha`, `lambda`, `mean_g`
and `sd_g` exceeds 1.05. The proportion of iterations at the maximum tree
depth and the cost are reported beside those of the settings of 0.3.0.

## Seeds

The stored fits have their own (20250807, 20250627). The fit with the
range from the true effects: 20250627. Equal standard errors: data
20261101, fit 20261102. Panels: data 20261009 + panel, fit 20261109 +
panel. ebnm is deterministic for the methods used here.

## Changes

1. Efron's estimator (2026-10-09, before any result of the comparison was
   seen). `deconv()` called with the estimates and `family = "Normal"`
   first counts them in 39 intervals. It is called with `P`, `Q` and
   `y = 1` instead: `P` is the normal likelihood of every estimate at
   every grid point and `Q` is the basis that `deconv()` builds for the
   normal family (`splines::ns(tau, M)`, its columns centered and scaled to
   length one), with `c0 = 1`. Both methods then use the same likelihood.
   `family = "Binomial"` is given with them: it selects the path of
   `deconv()` that takes one row of `P` for each observation, and has no
   other effect when `P` is given.
   The fit of bayesEfron for these data is
   `bayes_efron_fit(z, rep(1, 1500), M = 8, seed = 20261102, thin = 15)`,
   with 8 spline functions as for the simulated sites.
2. Fits of bayesEfron made for the comparison store 1 draw in 15, as the
   stored fits do (2026-10-09, before any result was seen).
3. Equal standard errors, repeated (2026-10-09, added after the result for
   the one data set was seen). In that data set the 90% intervals of
   bayesEfron contain 87.0% of the true effects and those of Efron's
   estimator 90.5%. One data set does not show whether that is a property
   of the method. `data-raw/benchmark/equal-se-run.R` repeats the
   comparison on the true effects of the 20 independent panels of 1,500
   sites, each with a new error of standard deviation 1 (data seed
   20261200 + panel, fit seed 20261300 + panel), with 6 and with 8 spline
   functions. Reported: the proportion of the true effects inside the 90%
   intervals, the mean length of those intervals and the root mean squared
   error for bayesEfron, its g plugged in and Efron's estimator, as means
   over the panels with their Monte Carlo standard errors. This is an
   addition; nothing fixed above is changed.
