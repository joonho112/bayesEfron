# Fitted Efron log-spline model objects

Objects of class `c("bef_fit_re", "bef_fit")`, as returned by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
contain the posterior draws, the data, the settings of the fit and
summaries of the draws (see Details). There are methods for printing and
summarizing a fit, for posterior means, intervals and quantiles of the
site effects, and for extracting the posterior draws. For plots,
predictions, diagnostics and model comparison, see
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md),
[`predict.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
and
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md).

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
coef(object, type = c("mean", "mode", "map"), ...)

# S3 method for class 'bef_fit_re'
confint(
  object,
  parm = NULL,
  level = 0.9,
  type = c("theta", "g"),
  method = c("exact", "draws"),
  interval = c("equal-tailed", "hpd"),
  ...
)

# S3 method for class 'bef_fit_re'
vcov(object, full = FALSE, ...)

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

# S3 method for class 'bef_fit_re'
quantile(
  x,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  reading = c("interpolated", "step"),
  ...
)

# S3 method for class 'summary.bef_fit'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_fit_re'
residuals(object, ...)

# S3 method for class 'bef_fit_re'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'bef_fit'
ndraws(x)

# S3 method for class 'bef_fit'
variables(x, ...)

# S3 method for class 'bef_fit'
AIC(object, ..., k = 2)

# S3 method for class 'bef_fit'
BIC(object, ...)
```

## Arguments

- x, object:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).
  For the `summary.bef_fit` methods, the result of
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Arguments passed on to
  [`format()`](https://rdrr.io/r/base/format.html) by the
  [`print()`](https://rdrr.io/r/base/print.html) methods.
  [`summary()`](https://rdrr.io/r/base/summary.html) throws an error if
  any are supplied; the other methods ignore them.

- level:

  Posterior probability of the intervals, a number between 0 and 1.

- type:

  For [`coef()`](https://rdrr.io/r/stats/coef.html), the summary to
  return: `"mean"` (the default) for the posterior means or `"mode"` for
  the posterior modes (`"map"`, the earlier name of `"mode"`, is also
  accepted). For [`confint()`](https://rdrr.io/r/stats/confint.html),
  the quantities for which intervals are returned: `"theta"` (the
  default) for the site effects or `"g"` for the mean, variance and
  standard deviation of \\g\\.

- parm:

  The sites or parameters for which
  [`confint()`](https://rdrr.io/r/stats/confint.html) returns intervals.
  For `type = "theta"`, a numeric vector of site positions in the data
  or a character vector of site labels. For `type = "g"`, one or more of
  `"mean_g"`, `"var_g"` and `"sd_g"`, or their positions. If `NULL` (the
  default), all are returned.

- method:

  How [`confint()`](https://rdrr.io/r/stats/confint.html) computes the
  limits for `type = "theta"`. `"exact"` (the default) takes them from
  the posterior distribution of the site effect, which averages over the
  draws of \\g\\. `"draws"` takes them as sample quantiles of the
  `theta_rep` draws and requires a fit made with `site_draws = TRUE`.

- interval:

  Type of interval returned by
  [`confint()`](https://rdrr.io/r/stats/confint.html). `"equal-tailed"`
  (the default) has its limits at the quantiles `(1 - level) / 2` and
  `1 - (1 - level) / 2`. `"hpd"` is the shortest interval with posterior
  probability `level`; if the distribution has one mode, it is the
  highest posterior density interval. For `method = "draws"` and for
  `type = "g"`, `"hpd"` gives the shortest interval that contains at
  least that proportion of the draws. The two types differ whenever the
  posterior distribution is not symmetric, and most when it has two
  modes.

- full:

  If `TRUE`, [`vcov()`](https://rdrr.io/r/stats/vcov.html) also returns
  the covariances between sites. If `FALSE` (the default), the matrix is
  diagonal.

- row.names:

  Optional character vector of row names for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), one
  per site.

- optional:

  Not used; it is an argument of the generic.

- use_cli:

  If `TRUE` and the cli package is installed, the first line of the
  output is bold and, in a printed summary, the section titles are
  colored. If `FALSE`, the output is plain text. If `NULL` (the
  default), the value of `getOption("bayesEfron.use_cli", TRUE)` is
  used.

- probs:

  Numeric vector of probabilities for
  [`quantile()`](https://rdrr.io/r/stats/quantile.html), each greater
  than 0 and less than 1.

- reading:

  How [`quantile()`](https://rdrr.io/r/stats/quantile.html) computes
  quantiles from point masses on the grid. `"interpolated"` (the
  default) spreads the probability of each grid point evenly over one
  grid spacing around it, and `"step"` leaves it on the grid point. See
  the Details of
  [`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md).

- nsim:

  Number of joint draws for
  [`simulate()`](https://rdrr.io/r/stats/simulate.html), a positive
  integer.

- seed:

  Integer seed for the draws of
  [`simulate()`](https://rdrr.io/r/stats/simulate.html), or `NULL` (the
  default) for none. If a seed is given, the state of R's random number
  generator is restored afterward.

- k:

  Not used.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns a list of
class `summary.bef_fit`. [`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return numeric
vectors with one element per site.
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
data frames, and [`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`quantile()`](https://rdrr.io/r/stats/quantile.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) return matrices.
The [`print()`](https://rdrr.io/r/base/print.html) methods return `x`
invisibly. The other methods are described in the two lists below:
"Methods (by generic)" for class `bef_fit`, and "Functions" for class
`bef_fit_re` and for a summary.

## Details

A fitted object is a list with the components `draws`, `posterior` and
`metadata`, and with `cmdstan_fit` if the model was fitted with
`keep_cmdstan_fit = TRUE`.

`draws` is a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
with the posterior draws of the following variables.

- `lp__`: The log density of the posterior distribution up to a
  constant, as reported by CmdStan.

- `alpha`: The spline coefficients.

- `lambda`: The precision of the prior for the spline coefficients.

- `log_g`: The log probabilities of \\g\\ at the grid points.

- `mean_g`, `var_g`, `sd_g`: The mean, variance and standard deviation
  of \\g\\.

- `effective_params`: The effective number of parameters defined below.

- `log_marginal_likelihood`: The marginal log likelihood of the
  estimates, with the site effects integrated out.

A fit made with `site_draws = TRUE` stores five more variables per site.

- `theta_mean`, `theta_sd`, `theta_map`: The mean, standard deviation
  and mode of the posterior distribution of \\\theta_i\\ given \\g\\.

- `theta_rep`: A draw from that distribution.

- `log_lik`: The log likelihood of the site's estimate.

Otherwise these variables are not stored, and the methods compute the
site summaries from the draws of \\g\\.

`posterior` is a list with the draws of `mean_g`, `var_g`, `sd_g`,
`effective_params` and `log_marginal_likelihood` as numeric vectors.

`metadata` is a list with the data, the settings of the fit and
summaries of the draws. It has the following components.

- `model_family`: The model family, currently always `"RE"`.

- `grid_method`: The method that chose the limits of the grid, or
  `"support"` if the limits were specified directly.

- `seed`: The seed used by the sampler.

- `cmdstan_version`: The version of CmdStan.

- `stan_file_sha256`: The SHA-256 hash of the Stan file.

- `data_list`: The data passed to Stan, a list with the estimates
  `theta_hat`, the standard errors `sigma`, the grid points `grid`, the
  spline basis `B`, and the numbers of sites, grid points and spline
  functions (`K`, `L` and `M`).

- `runtime_seconds`: The time taken by sampling, in seconds.

- `theta_summary`: The table of site effects, a data frame with one row
  per site. Its columns are described below.

- `mean_g_summary`, `var_g_summary`, `effective_params_summary`,
  `log_marginal_likelihood_summary`: Lists with the mean, the standard
  deviation and the 5%, 50% and 95% quantiles of the draws.

More information is stored in attributes of `metadata` and can be
extracted with [`attr()`](https://rdrr.io/r/base/attr.html), e.g.,
`attr(fit$metadata, "grid")`.

- `sampler_settings`: A list with the values of `chains`,
  `parallel_chains`, `iter_warmup`, `iter_sampling`, `adapt_delta`,
  `max_treedepth`, `seed`, `refresh` and `init` passed to CmdStan.

- `grid`: The grid of the fit, as returned by
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md).
  It can be supplied as `grid` in another fit.

- `storage`: A list recording whether the site-level draws are stored,
  the number of chains, the number of iterations stored per chain, and
  the value of `thin`.

- `package_version`, `cmdstanr_version`: The versions of bayesEfron and
  cmdstanr used for the fit.

`theta_summary`, the table of site effects, has the columns listed
below. The columns from `mean` to `mode` summarize the posterior
distribution of \\\theta_i\\, which averages over the posterior draws of
\\g\\ (Lee and Sui, 2025, Section 4.3).

- `site`: The site number, in the order of the data.

- `label`: The site label, or the site number if the data had no labels.

- `estimate`, `se`: The estimate and its standard error, as supplied.

- `mean`: The posterior mean.

- `sd`: The posterior standard deviation. The posterior variance is the
  variance of \\\theta_i\\ given \\g\\, averaged over the draws of
  \\g\\, plus the variance over those draws of the mean of \\\theta_i\\
  given \\g\\, so `sd` includes the uncertainty about \\g\\.

- `lower`, `upper`: The 5% and 95% quantiles of the posterior
  distribution, i.e., the limits of an equal-tailed 90% interval. They
  are quantiles of the distribution itself and not sample quantiles of
  site-level draws. For a fit with point masses they are computed with
  the `"interpolated"` reading (described under `reading`), so they
  generally fall between grid points.

- `mode`: The posterior mode. For point masses this is the grid point
  with the largest posterior probability, and for normal components the
  value at which the posterior density is largest.

- `hpdi_lower`, `hpdi_upper`, `map`: Deprecated columns that repeat
  `lower`, `upper` and `mode`. They will be removed.

Like all site summaries, the posterior summaries in this table are
computed from the stored draws of \\g\\ and are subject to the Monte
Carlo error of those draws.
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
returns the same table for any level, with additional columns.

`effective_params` is the sum over sites of the ratio of the posterior
variance of \\\theta_i\\ given \\g\\ to \\\sigma_i^2\\. In the normal
random-effects model this ratio is the weight given to the site's own
estimate, so the sum is small when the estimates are shrunk strongly and
close to the number of sites when there is little shrinkage. Note that
`p_loo` in
[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md),
which the loo package also calls the effective number of parameters, is
a different quantity.

In the printed summary, the mean, variance and standard deviation of
\\g\\ are labeled "Prior g", because \\g\\ is the prior distribution of
each \\\theta_i\\.

A fit saved by an earlier version of the package may contain a table of
site effects or a value of `effective_params` computed under other
definitions. The methods then recompute them from the stored draws of
\\g\\ and issue a warning that names the recomputed quantities. The
saved object is not changed.

## Methods (by generic)

- `print(bef_fit)`: Prints the number of sites, the number of draws, the
  worst value of each sampler diagnostic, the range of the standard
  errors, the grid and the method that chose its range, and the
  functions to call next.

- `summary(bef_fit)`: Returns a list of class `summary.bef_fit` with the
  components `prior_summary` (the posterior means of the mean, variance
  and standard deviation of \\g\\) and `diagnostics` (the sampler
  diagnostics; see
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)).

- `nobs(bef_fit)`: Returns the number of sites.

- `logLik(bef_fit)`: Returns the posterior mean of
  `log_marginal_likelihood`, the marginal log likelihood with the site
  effects integrated out, as a `logLik` object. Its degrees of freedom
  are the posterior mean of `effective_params`. Note that this is not a
  maximized log likelihood, so
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) throw an error. Use
  [`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)
  to compare fits.

- `as_draws(bef_fit)`: Returns the posterior draws as a
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html).

- `format(bef_fit)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows as a character
  vector.

- `ndraws(bef_fit)`: Returns the number of posterior draws stored in the
  fit, summed over the chains.

- `variables(bef_fit)`: Returns the names of the variables in the draws.

- `AIC(bef_fit)`: [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) throw an error, for the
  reason given under [`logLik()`](https://rdrr.io/r/stats/logLik.html).

## Functions

- `summary(bef_fit_re)`: Adds the component `theta_summary` to the
  summary: the table of site effects described in Details, with
  equal-tailed intervals at `level`.

- `print(summary.bef_fit)`: Prints the summary.

- `coef(bef_fit_re)`: Returns the posterior means or the posterior modes
  of the site effects (the column `mean` or `mode` of `theta_summary`)
  as a numeric vector. Its names are the site labels, or the site
  numbers if the data had no labels.

- `confint(bef_fit_re)`: Returns a data frame of posterior intervals
  with the attribute `level`. For `type = "theta"` it has one row per
  site and the columns `site`, `mean`, `lower` and `upper`. For
  `type = "g"` it has one row each for the mean, variance and standard
  deviation of \\g\\ and the columns `parameter`, `mean`, `lower` and
  `upper`. The limits are computed as specified by `method` and
  `interval`. The column `point` and, for `type = "g"`, the column
  `site` repeat `mean` and `parameter`. They are deprecated and will be
  removed.

- `vcov(bef_fit_re)`: Returns the posterior variances of the site
  effects (the squares of the column `sd` of `theta_summary`) as a
  diagonal matrix. If `full = TRUE`, the matrix also contains the
  posterior covariances between sites. The site effects are independent
  given \\g\\ but not in the posterior distribution, because all of them
  depend on \\g\\. Their covariances are the covariances, over the draws
  of \\g\\, of the means of the effects given \\g\\. The diagonal matrix
  is sufficient when sites are considered one at a time. Use
  `full = TRUE` for the variance of a sum or a difference of site
  effects.

- `fitted(bef_fit_re)`: Returns the posterior means of the site effects,
  as `coef(object)` does.

- `as.data.frame(bef_fit_re)`: Returns `theta_summary`, the table of
  site effects described in Details, with 90% intervals.

- `quantile(bef_fit_re)`: Returns quantiles of the posterior
  distributions of the site effects as a matrix with one row per site
  and one column per probability in `probs`.

- `format(summary.bef_fit)`: Returns the lines of the printed summary as
  a character vector.

- `residuals(bef_fit_re)`: Returns the estimates minus the posterior
  means of the site effects, as a numeric vector named by the site
  labels.

- `simulate(bef_fit_re)`: Draws the site effects jointly from their
  posterior distribution. For each of the `nsim` draws, one posterior
  draw of \\g\\ is taken and every site effect is drawn from its
  distribution given that \\g\\. Returns a matrix with one row per draw
  and one column per site. Each row is one joint draw, so functions of
  several sites (e.g., ranks or the largest effect) can be computed row
  by row.

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
fit <- krw_fit
fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 97
#> Grid method: pilot
#> Runtime: 0.7197 sec
#> Draws: 800 (4 chains, one iteration in 15 kept); site-level draws not kept
#> Diagnostics: Rhat 1.002; ESS bulk 1731; ESS tail 1780; divergences 0; treedepth hits 0
#> Standard errors: 0.005344 to 0.0392 (largest / smallest 7.335)
#> Grid: 101 points from -0.0829 to 0.1585; smallest standard error 2.213 grid spacings
#> Grid range: rule "pilot", not narrowed (the range of "paper_realdata"); point masses; 6 spline functions
#> Next: summary() and site_summary() for the sites, g_summary() for the distribution of effects, plot().
summary(fit)$prior_summary
#> $mean
#> [1] 0.01660178
#> 
#> $var
#> [1] 0.0003076948
#> 
#> $sd
#> [1] 0.01716216
#> 

# Posterior means, next to the estimates they are shrunk from.
head(cbind(estimate = krw_discrimination$yi, posterior_mean = coef(fit)))
#>                        estimate posterior_mean
#> Rite Aid            0.046957123    0.030481747
#> KFC                 0.022000000    0.017864935
#> Victoria's Secret   0.042161613    0.023347817
#> United Rentals      0.005708306    0.010719906
#> Cardinal Health     0.034077110    0.021096229
#> Charter / Spectrum -0.010182468    0.001675308

# Intervals for three sites, and for the mean and spread of g.
confint(fit, parm = c(29, 11, 65), level = 0.95)
#> 95% posterior intervals (equal-tailed)
#>                   site         mean        lower       upper        point
#> Quest Diagnostics   29 -0.004794688 -0.014809407 0.004878364 -0.004794688
#> J.C. Penney         11  0.021602438  0.001058169 0.045289468  0.021602438
#> Walmart             65  0.023868227 -0.003405634 0.065165090  0.023868227
confint(fit, type = "g")
#> 90% posterior intervals (equal-tailed)
#>   parameter         mean        lower        upper   site        point
#> 1    mean_g 0.0166017761 0.0130811925 0.0206491286 mean_g 0.0166017761
#> 2     var_g 0.0003076948 0.0001399332 0.0005738865  var_g 0.0003076948
#> 3      sd_g 0.0171621568 0.0118293270 0.0239559271   sd_g 0.0171621568

# Posterior standard deviations and posterior modes, with the standard
# errors of the estimates.
head(cbind(se = krw_discrimination$sei, as.data.frame(fit)[c("mean", "sd", "mode")]))
#>           se        mean          sd        mode
#> 1 0.01619358 0.030481747 0.012967949 0.028158038
#> 2 0.01530472 0.017864935 0.010525940 0.016086924
#> 3 0.02296031 0.023347817 0.013848996 0.020915369
#> 4 0.01504751 0.010719906 0.009844949 0.008844255
#> 5 0.02145421 0.021096229 0.012896169 0.018501147
#> 6 0.01209059 0.001675308 0.008784082 0.001601586

# Posterior variances as a diagonal matrix.
diag(vcov(fit))[1:3]
#>          Rite Aid               KFC Victoria's Secret 
#>      0.0001681677      0.0001107954      0.0001917947 
```
