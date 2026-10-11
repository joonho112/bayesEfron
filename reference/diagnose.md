# MCMC convergence and grid diagnostics for a fitted model

`diagnose()` returns the MCMC diagnostics of a fit (R-hat, effective
sample sizes, divergent transitions, tree depth and E-BFMI) and flags
those that indicate a problem. The diagnostics were computed by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
when the model was fitted, so no sampling is done. `diagnose()` also
checks the grid, in particular whether the result for any site depends
on where the grid ends.

## Usage

``` r
diagnose(fit, ...)

# Default S3 method
diagnose(fit, ...)

# S3 method for class 'bef_fit'
diagnose(fit, ...)

# S3 method for class 'bef_diagnostic'
print(x, ...)

# S3 method for class 'bef_diagnostic'
summary(object, ...)

# S3 method for class 'bef_diagnostic'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'bef_diagnostic'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_fit'
nuts_params(object, pars = NULL, ...)

# S3 method for class 'bef_fit'
log_posterior(object, ...)

# S3 method for class 'bef_fit'
rhat(object, pars = NULL, ...)

# S3 method for class 'bef_fit'
neff_ratio(object, pars = NULL, ...)
```

## Arguments

- fit:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- ...:

  For `diagnose()`, these dots are for future extensions and must be
  empty. The [`print()`](https://rdrr.io/r/base/print.html) method
  passes them on to [`format()`](https://rdrr.io/r/base/format.html);
  [`format()`](https://rdrr.io/r/base/format.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) ignore them.

- x, object:

  A `bef_diagnostic` object returned by `diagnose()`. For the bayesplot
  methods (`nuts_params()`, `log_posterior()`, `rhat()` and
  `neff_ratio()`), `object` is a `bef_fit` object.

- row.names, optional:

  Ignored. They are arguments of the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic.

- use_cli:

  If `TRUE` and the cli package is installed, the first line is
  formatted in bold; `FALSE` gives plain text. If `NULL` (the default),
  the value of `getOption("bayesEfron.use_cli", TRUE)` is used.

- pars:

  An optional character vector. For `nuts_params()`, the names of the
  sampler quantities to return; for `rhat()` and `neff_ratio()`, the
  names of the variables. If `NULL` (the default), all are returned.

## Value

An object of class `bef_diagnostic`, which is a list with the following
components.

- `rhat`:

  R-hat for each variable in the draws, a named numeric vector.

- `ess_bulk`, `ess_tail`:

  The bulk and tail effective sample sizes, in the same form.

- `divergences`, `max_treedepth`:

  The number of divergent transitions and the number of iterations that
  reached the maximum tree depth, each summed over the chains.

- `ebfmi`:

  The E-BFMI of each chain. Like the two counts, it is `NA` if CmdStan
  did not report it.

- `sampler_diagnostics_failed`:

  Names of the diagnostics that fail a check.

- `sampler_diagnostics_warned`:

  Names of the diagnostics that lead to a warning. `"grid_ends"` is
  among them if an end of the grid has the status `"warning"`.

- `diagnostic_skipped`:

  Names of the diagnostics that could not be computed or were not
  checked. These three character vectors are empty if the run has at
  least 400 draws and there are no problems.

- `grid`:

  A list with the grid diagnostics (see below), or `NULL` if they cannot
  be computed for the fit.

- `model_family`, `stan_file_sha256`, `runtime_seconds`,
  `effective_params_summary`:

  Values copied from `fit$metadata` (see
  [bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)).

The list `grid` has the following components.

- `ends`:

  A data frame with one row per end of the grid (see below).

- `share`:

  The outer fraction of the grid's range that is examined at each end
  (0.02).

- `threshold`:

  The posterior probability above which an end is flagged (0.05).

- `sd_g`:

  The posterior mean and the 5% and 95% quantiles of the standard
  deviation of \\g\\.

- `sd_g_prior`:

  The 5%, 50% and 95% quantiles of the prior distribution of that
  standard deviation on this grid.

- `occupied`:

  The proportion of the grid's range that contains the central 99.9% of
  the estimate of \\g\\.

- `tau_c_over_sd_g`:

  The standard deviation of a normal component divided by the standard
  deviation of \\g\\ (`NA` for point masses).

The data frame `grid$ends` has the following columns.

- `end`:

  `"lower"` or `"upper"`.

- `value`:

  The grid point at that end.

- `set_by`:

  How the end was set: `"support"`, `"estimates"`, `"pilot"` or
  `"true effects"`.

- `mass`:

  The largest posterior probability that any site has in the outer 2% of
  the range at that end.

- `site`, `label`:

  The position in the data and the label of that site.

- `estimate_beyond`:

  How far the estimate of that site lies beyond the end, in units of its
  standard error (negative if it lies inside the grid).

- `sites_over`:

  The number of sites for which this probability exceeds 0.05.

- `limits_at_end`:

  The number of sites whose 90% interval has a limit at the grid point
  at that end when `reading = "step"` is used, because that grid point
  alone has at least the tail probability of the interval.

- `status`:

  `"ok"`, `"warning"` or `"note"`.

## Details

A diagnostic fails its check and is named in
`sampler_diagnostics_failed` in two cases.

- An R-hat value exceeds 1.05.

- There is a divergent transition.

A diagnostic leads to a warning and is named in
`sampler_diagnostics_warned` in the following cases.

- An R-hat value is between 1.01 and 1.05.

- A bulk or tail effective sample size is below 400.

- An iteration reached the maximum tree depth.

- The energy Bayesian fraction of missing information (E-BFMI) of a
  chain is below 0.3.

R-hat below 1.01 and effective sample sizes above 400 are the
recommendations of Vehtari et al. (2021), and 0.3 is the value below
which CmdStan warns about E-BFMI. The other thresholds are choices of
this package. See
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
for what to do about each kind of warning.

The sampler diagnostics describe the whole run, even if the fit was
thinned. If the run has fewer than 400 draws in total, R-hat and the
effective sample sizes are too unreliable to act on. Their values are
still returned, but they are not compared with the thresholds, and
`diagnostic_skipped` contains `"rhat_check"`, `"ess_bulk_check"` and
`"ess_tail_check"`. Divergent transitions, the maximum tree depth and
E-BFMI are checked regardless of the length of the run.

R-hat and the effective sample sizes are `NA` for a quantity that has
the same value in every draw, e.g., the posterior mode of a site that
never leaves one grid point in a fit with `site_draws = TRUE`. If a
quantity is constant within each chain but differs between chains, the
chains have not mixed and its R-hat is infinite. That R-hat is stored as
`NA`, and `"rhat"` is named in `sampler_diagnostics_failed` even for a
run of fewer than 400 draws. If R-hat or an effective sample size cannot
be computed at all, or an effective sample size is infinite, all values
of that diagnostic are `NA` and it is named in `diagnostic_skipped`.

The `grid` component contains the grid diagnostics. The main one is a
check of the ends of the grid, i.e., its lower and upper limits. With
point masses, the model places no effect beyond them. For each end,
`grid$ends` gives the largest posterior probability that any site has in
the outer 2% of the grid's range (`mass`). If this probability exceeds
0.05, the result for that site depends on where the grid ends, and the
`status` of the end depends on how the end was set (`set_by`). Otherwise
the status is `"ok"`.

- The status is `"warning"` if the end was specified in `support` or
  follows from the estimates alone (`set_by` is `"support"` or
  `"estimates"`).
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  then issues a warning of class `bef_warn_grid_ends` that names the
  site, and `"grid_ends"` is added to `sampler_diagnostics_warned`,
  although it concerns the grid and not the sampler. A wider range may
  change the result, and
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  shows whether it does.

- The status is `"note"` if the `"pilot"` method moved the end in
  (`set_by` is `"pilot"`). The pilot fit, which used a wider range,
  placed at most 0.05% of \\g\\ beyond that end. A site with a large
  probability there has an estimate beyond the end of \\g\\ and is
  shrunk to it, so the probability is reported without a warning. The
  same applies to a range around the true effects of a simulation
  (`set_by` is `"true effects"`).

See the examples of
[`plot.bef_grid_comparison()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid_comparison.md)
for a comparison in which the narrowest range gives a probability above
0.05.

No single number shows whether a range is wider than necessary. If there
are few sites, the data carry little information about the parts of the
range that contain no estimates, and \\g\\ is determined by the prior
there. The `grid` component therefore gives the posterior mean and 90%
interval of the standard deviation of \\g\\ (`sd_g`) next to the 5%, 50%
and 95% quantiles of its prior distribution on the same grid
(`sd_g_prior`). If the two are close, the data have hardly determined
the tails of \\g\\, and it is worth comparing fits with other ranges.
The proportion of the grid's range that contains the central 99.9% of
the estimate of \\g\\ (`occupied`) is informative only if there are many
sites.

## Functions

- `print(bef_diagnostic)`: Prints the worst value of each diagnostic.

- `summary(bef_diagnostic)`: Returns a list in which `rhat`, `ess_bulk`
  and `ess_tail` are reduced to their worst value (the largest R-hat and
  the smallest effective sample sizes), with the name of the variable
  that has it.

- `as.data.frame(bef_diagnostic)`: Returns a data frame with one row per
  diagnostic and the columns `diagnostic`, `value` (the worst value),
  `threshold`, `rule` (whether `value` should be `"at most"` or
  `"at least"` the threshold) and `status` (`"ok"`, `"warning"`,
  `"failed"` or `"not checked"`). The last row, `grid_ends`, is for the
  end of the grid at which a site has the largest posterior probability
  in the outer 2% of the range. Its status is `"note"` if that
  probability exceeds the threshold at an end set by the `"pilot"`
  method or from the true effects.

- `format(bef_diagnostic)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows as a character
  vector.

- `nuts_params(bef_fit)`: Methods for bayesplot generics, so that the
  bayesplot plots of sampler diagnostics can be made from a fit.
  `nuts_params()` returns the divergence indicator, the tree depth, the
  energy, the acceptance statistic, the step size and the number of
  leapfrog steps for every stored iteration, and `log_posterior()`
  returns the log density for every stored iteration. `rhat()` and
  `neff_ratio()` return, for each variable, R-hat and the bulk effective
  sample size divided by the number of iterations of the whole run. In a
  thinned fit (see `thin` in
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md))
  the stored iterations are not consecutive, so a plot that compares
  each iteration with the previous one, such as
  [`bayesplot::mcmc_nuts_energy()`](https://mc-stan.org/bayesplot/reference/MCMC-nuts.html),
  should not be made from it. `nuts_params()` throws an error if the fit
  does not contain the sampler quantities.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis*, 16(2).
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md),
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)

## Examples

``` r
diag <- diagnose(krw_fit)
diag
#> <bef_diagnostic>
#> Model family: RE
#> Rhat max: 1.002 (log_g[97])
#> ESS bulk min: 1731 (lp__)
#> ESS tail min: 1780 (log_g[2])
#> Divergences: 0
#> Max treedepth hits: 0
#> E-BFMI min: 0.5682
#> Runtime: 0.7197 sec
#> Grid:
#>   lower end -0.0829: at most 0% of a site's posterior in the outer 2% (site "Safeway")
#>   upper end 0.1585: at most 0.1% of a site's posterior in the outer 2% (site "Genuine Parts (Napa Auto)")
#>   sd of g: posterior 0.01716 (0.01183 to 0.02396); prior 0.06991 (0.05814 to 0.07644)
#>   the central 99.9% of g spans 98% of the grid
#> Effective params: mean 46.67; sd 6.378

# The largest R-hat and the variable that has it.
summary(diag)$rhat
#> $value
#> [1] 1.001951
#> 
#> $index
#> [1] 105
#> 
#> $variable
#> [1] "log_g[97]"
#> 

# The three variables with the smallest bulk effective sample size.
head(sort(diag$ess_bulk), 3)
#>     lp__   lambda log_g[6] 
#> 1731.204 2082.036 2384.530 
```
