# Fitted Efron log-spline model objects

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
returns an object of class `bef_fit_re`, which inherits from `bef_fit`.
This page describes what the object holds and the methods that summarize
it.
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md),
[`predict.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
and
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
have pages of their own.

## Usage

``` r
# S3 method for class 'bef_fit'
print(x, ...)

# S3 method for class 'bef_fit'
summary(object, level = 0.9, ...)

# S3 method for class 'bef_fit_re'
summary(object, level = 0.9, ...)

# S3 method for class 'summary.bef_fit'
print(x, ...)

# S3 method for class 'bef_fit_re'
coef(object, type = c("mean", "map"), ...)

# S3 method for class 'bef_fit_re'
confint(object, parm = NULL, level = 0.9, type = c("theta", "g"), ...)

# S3 method for class 'bef_fit_re'
vcov(object, ...)

# S3 method for class 'bef_fit_re'
fitted(object, ...)

# S3 method for class 'bef_fit_re'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'bef_fit'
nobs(object, ...)

# S3 method for class 'bef_fit'
logLik(object, ...)

# S3 method for class 'bef_fit'
as_draws(x, ...)

# S3 method for class 'bef_fit'
format(x, ..., use_cli = NULL)

# S3 method for class 'summary.bef_fit'
format(x, ..., use_cli = NULL)
```

## Arguments

- x, object:

  A `bef_fit` object from
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
  or for the `summary.bef_fit` methods the result of
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) methods,
  arguments passed on to
  [`format()`](https://rdrr.io/r/base/format.html). The other methods do
  not use them.

- level:

  Number between 0 and 1, the probability content of the intervals.
  Defaults to `0.9`.

- type:

  For [`coef()`](https://rdrr.io/r/stats/coef.html), `"mean"` (the
  default) for the posterior means or `"map"` for the posterior modes.
  For [`confint()`](https://rdrr.io/r/stats/confint.html), `"theta"`
  (the default) for intervals for the site effects or `"g"` for
  intervals for the mean, variance and standard deviation of \\g\\.

- parm:

  The rows to return: site numbers for `type = "theta"`, and for
  `type = "g"` one or more of `"mean_g"`, `"var_g"` and `"sd_g"` or
  their positions. `NULL` (the default) returns all.

- row.names:

  Character vector of row names, one for each site, or `NULL` (the
  default).

- optional:

  Not used; it is an argument of the generic.

- use_cli:

  Logical or `NULL`. With `TRUE` and the cli package installed, the
  first line is set in bold and section titles in color; `FALSE` gives
  plain text. `NULL` (the default) takes the value of
  `getOption("bayesEfron.use_cli", TRUE)`.

## Details

The object is a list with three elements.

`draws` is a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
with the posterior draws of every quantity in the Stan model: the spline
coefficients `alpha` and the precision `lambda`; `log_g`, the log
probabilities of \\g\\ on the grid; `mean_g`, `var_g` and `sd_g`, the
mean, variance and standard deviation of \\g\\; for each site the mean
`theta_mean`, standard deviation `theta_sd` and mode `theta_map` of the
posterior distribution of \\\theta_i\\ given \\g\\, a draw `theta_rep`
from that distribution and the log likelihood `log_lik`;
`effective_params`; and `log_marginal_likelihood`, the sum of `log_lik`
over the sites.

`posterior` is a list with the draws of `mean_g`, `var_g`, `sd_g`,
`effective_params` and `log_marginal_likelihood` as numeric vectors.

`metadata` is a list with the settings of the fit and summaries of the
draws:

- `model_family` and `grid_method`, as given to
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- `seed`, the seed that the sampler used; `cmdstan_version`; and
  `stan_file_sha256`, a hash of the Stan file.

- `data_list`, the data passed to Stan: the estimates `theta_hat`, the
  standard errors `sigma`, the grid points `grid`, the spline basis `B`,
  and the numbers of sites, grid points and spline functions, `K`, `L`
  and `M`.

- `runtime_seconds`, the time that sampling took.

- `theta_summary`, a data frame with one row for each site, described
  below.

- `mean_g_summary`, `var_g_summary`, `effective_params_summary` and
  `log_marginal_likelihood_summary`, each a list with the mean, the
  standard deviation and the 5%, 50% and 95% quantiles of the draws.

The settings of the sampler are attributes of `metadata`.
`attr(fit$metadata, "sampler_settings")` is a list with the values of
`chains`, `parallel_chains`, `iter_warmup`, `iter_sampling`,
`adapt_delta`, `max_treedepth`, `seed`, `refresh` and `init` that were
passed to CmdStan, and the attributes `package_version` and
`cmdstanr_version` give the versions of bayesEfron and cmdstanr that
made the fit. A fit saved by an earlier version of the package does not
have these attributes.

If the model was fitted with `keep_cmdstan_fit = TRUE` there is a fourth
element, `cmdstan_fit`.

In `theta_summary`, `site` numbers the sites in the order of the data.
The other columns summarize the posterior distribution of \\\theta_i\\,
which averages over the posterior draws of \\g\\ (Lee and Sui, 2025,
Section 4.3). `mean` is its mean and `sd` its standard deviation. The
posterior variance is the variance of \\\theta_i\\ given \\g\\, averaged
over the draws of \\g\\, plus the variance over those draws of the mean
of \\\theta_i\\ given \\g\\, so that `sd` includes the uncertainty about
\\g\\. `hpdi_lower` and `hpdi_upper` are the 5% and 95% quantiles of the
posterior draws `theta_rep`; in spite of their names they are the limits
of an equal-tailed interval and not of a highest posterior density
interval. The draws take values on the grid, but a sample quantile is
interpolated between two adjacent values of the ordered draws, so a
limit can lie between two grid points. `map` is the posterior mode, the
grid point to which the posterior distribution gives the largest
probability.

`effective_params` is the sum over the sites of the ratio of the
posterior variance of \\\theta_i\\ given \\g\\ to \\\sigma_i^2\\. In the
normal random-effects model that ratio is the weight given to the site's
own estimate, so the sum is small when the estimates are shrunk strongly
and near the number of sites when they are hardly shrunk. It is a
different quantity from the `p_loo` of
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md),
which the loo package also calls an effective number of parameters.

\\g\\ is the prior distribution of each \\\theta_i\\, and the printed
summary labels its mean, variance and standard deviation "Prior g".

The definitions of `sd`, `map` and `effective_params` have changed
between versions of the package, and a fit records which definitions its
summaries follow, in the attribute `summary_definition_version` of
`metadata`. For a fit that was saved by an earlier version, the methods
on this page compute `sd` and `map` again from the stored draws, and
`effective_params` too if it was stored under its first definition. They
give a warning that names the quantities that changed and leave the
saved object as it is. The result of
[`summary()`](https://rdrr.io/r/base/summary.html) saved by an earlier
version holds no draws and cannot be brought up to date;
[`print()`](https://rdrr.io/r/base/print.html) stops and asks for
[`summary()`](https://rdrr.io/r/base/summary.html) to be called again on
the fit.
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`predict.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md)
and
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
need quantities that the first released version did not store, and stop
for a fit saved by it.

## Methods (by generic)

- `print(bef_fit)`: Prints the number of sites, the grid method, the
  running time and the worst value of each convergence diagnostic.

- `summary(bef_fit)`: Returns a list of class `summary.bef_fit` with
  `prior_summary`, the posterior means of the mean, variance and
  standard deviation of \\g\\, and `diagnostics`, the convergence
  diagnostics.

- `nobs(bef_fit)`: Returns the number of sites.

- `logLik(bef_fit)`: Returns the posterior mean of
  `log_marginal_likelihood`, the log likelihood with the site effects
  summed out, as a `logLik` object whose degrees of freedom are the
  posterior mean of `effective_params`. It is not a maximized log
  likelihood, so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) computed from it do not
  have their usual meaning; to compare fits, use
  [`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md).

- `as_draws(bef_fit)`: Returns the draws as a
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html).

- `format(bef_fit)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows, as a character
  vector.

## Functions

- `summary(bef_fit_re)`: Adds `theta_summary` to the summary, with the
  columns `hpdi_lower` and `hpdi_upper` computed for `level`.

- `print(summary.bef_fit)`: Prints the summary.

- `coef(bef_fit_re)`: Returns the posterior means or the posterior modes
  of the site effects, the column `mean` or `map` of `theta_summary`, as
  a numeric vector named by site number.

- `confint(bef_fit_re)`: Returns a data frame with the columns `site`,
  `lower`, `upper` and `point`. For `type = "theta"` the limits are
  equal-tailed quantiles of the `theta_rep` draws and `point` is the
  posterior mean. For `type = "g"` the rows are the mean, variance and
  standard deviation of \\g\\.

- `vcov(bef_fit_re)`: Returns a diagonal matrix with the posterior
  variances of the site effects, the squares of the column `sd` of
  `theta_summary`. The site effects are not independent in the posterior
  distribution, because all of them depend on \\g\\, and their
  covariances are not computed. The matrix therefore serves for one site
  at a time and not for the variance of a sum or a difference of site
  effects.

- `fitted(bef_fit_re)`: Returns the posterior means of the site effects,
  as `coef(object)` does.

- `as.data.frame(bef_fit_re)`: Returns `theta_summary`.

- `format(summary.bef_fit)`: Returns the lines of the printed summary.

## References

Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
2639. [doi:10.3390/math13162639](https://doi.org/10.3390/math13162639)

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md),
[`predict.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md),
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)

## Examples

``` r
fit <- raudenbush_fit
fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 19
#> Grid method: paper_realdata
#> Runtime: 1.351 sec
#> Diagnostics: Rhat 1.007; ESS bulk 676.6; ESS tail 888.8; divergences 0; treedepth hits 1
#> Warnings: max_treedepth
#> Use summary() for posterior summaries.
summary(fit)
#> <summary.bef_fit>
#> 
#> Prior g:
#>   mean: 0.1145
#>   var:  0.1094
#>   sd:   0.3009
#> 
#> Diagnostics:
#>   Rhat (max):        1.007
#>   ESS bulk (min):    676.6
#>   ESS tail (min):    888.8
#>   Divergences:       0
#>   Treedepth hits:    1
#>   Effective params:  mean 10.29, sd 2.84
#>   Log marginal lik.: mean -4.047, sd 1.439
#>   Runtime:           1.351 sec
#>   Warnings:          max_treedepth
#> 
#> Theta summary:
#>   site      mean        sd     lower     upper       map
#>   1      0.03805   0.09832     -0.11      0.19      0.04
#>   2       0.0874    0.1121     -0.08      0.28      0.07
#>   3     -0.03442    0.1188     -0.23      0.16     -0.02
#>   4        0.635    0.4979      0.04      1.48      0.19
#>   5       0.1063     0.186     -0.14    0.4015      0.07
#>   6     -0.02282   0.08604     -0.17       0.1     -0.02
#>   14      0.1073    0.1658     -0.11       0.4      0.07
#>   15    -0.05756    0.1177     -0.29    0.1015     -0.05
#>   16    0.001513    0.1169      -0.2      0.19      0.01
#>   17      0.1979    0.1221      0.01       0.4      0.19
#>   18       0.063   0.08101     -0.08      0.19      0.07
#>   19   -0.0006152    0.1197      -0.2      0.19      0.01
#>   ... 7 sites omitted ...

# Posterior means, next to the estimates they are shrunk from.
cbind(estimate = raudenbush1985$yi, posterior_mean = coef(fit))
#>    estimate posterior_mean
#> 1      0.03   0.0380499230
#> 2      0.12   0.0873985077
#> 3     -0.14  -0.0344193518
#> 4      1.18   0.6349620514
#> 5      0.26   0.1062798275
#> 6     -0.06  -0.0228176947
#> 7     -0.02   0.0031007250
#> 8     -0.32  -0.0749462778
#> 9      0.27   0.1627523821
#> 10     0.80   0.4652042240
#> 11     0.54   0.2175916906
#> 12     0.18   0.1011498959
#> 13    -0.02   0.0396075294
#> 14     0.23   0.1073244942
#> 15    -0.18  -0.0575610879
#> 16    -0.06   0.0015127205
#> 17     0.30   0.1979408996
#> 18     0.07   0.0630018949
#> 19    -0.07  -0.0006152146

# Intervals for three sites, and for the mean and spread of g.
confint(fit, parm = c(4, 10, 18), level = 0.95)
#>   site lower upper      point
#> 1    4 -0.02  1.66 0.63496205
#> 2   10  0.04  1.12 0.46520422
#> 3   18 -0.08  0.22 0.06300189
confint(fit, type = "g")
#>     site       lower     upper     point
#> 1 mean_g 0.003896369 0.2499901 0.1144608
#> 2  var_g 0.010904857 0.3041748 0.1093714
#> 3   sd_g 0.104426290 0.5515191 0.3008526

# Posterior standard deviations and posterior modes, with the standard
# errors of the estimates.
head(cbind(se = raudenbush1985$sei, as.data.frame(fit)[c("mean", "sd", "map")]))
#>          se        mean         sd   map
#> 1 0.1249000  0.03804992 0.09831588  0.04
#> 2 0.1469694  0.08739851 0.11205118  0.07
#> 3 0.1670329 -0.03441935 0.11878892 -0.02
#> 4 0.3729611  0.63496205 0.49794655  0.19
#> 5 0.3690528  0.10627983 0.18595126  0.07
#> 6 0.1029563 -0.02281769 0.08604151 -0.02

# Posterior variances as a diagonal matrix.
diag(vcov(fit))[1:3]
#>           1           2           3 
#> 0.009666013 0.012555467 0.014110806 
```
