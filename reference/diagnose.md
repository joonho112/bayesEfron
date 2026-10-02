# Convergence diagnostics of a fitted model

`diagnose()` collects the diagnostics of the sampler that
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
computed when the model was fitted and names those that fail a check. It
computes nothing new from the draws.

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
format(x, ..., use_cli = NULL)
```

## Arguments

- fit:

  A `bef_fit` object from
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- ...:

  For `diagnose()`, these dots are for future extensions and must be
  empty. The [`print()`](https://rdrr.io/r/base/print.html) method
  passes them on to [`format()`](https://rdrr.io/r/base/format.html);
  [`format()`](https://rdrr.io/r/base/format.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) do not use them.

- x, object:

  A `bef_diagnostic` object.

- use_cli:

  Logical or `NULL`. With `TRUE` and the cli package installed, the
  first line is set in bold and section titles in color; `FALSE` gives
  plain text. `NULL` (the default) takes the value of
  `getOption("bayesEfron.use_cli", TRUE)`.

## Value

A list of class `bef_diagnostic`. `rhat`, `ess_bulk` and `ess_tail` are
numeric vectors with one value for each quantity in the draws.
`divergences` and `max_treedepth` are the numbers of divergent
transitions and of iterations that reached the maximum tree depth,
totaled over the chains, and `ebfmi` has one value for each chain; the
three are `NA` if CmdStan did not report them. `model_family`,
`stan_file_sha256`, `runtime_seconds` and `effective_params_summary`
repeat values from `fit$metadata`. `sampler_diagnostics_failed`,
`sampler_diagnostics_warned` and `diagnostic_skipped` are character
vectors that name the diagnostics that fail a check, those that deserve
a warning, and those that could not be computed or were not checked (see
Details). All three are empty for a fit of 400 or more draws without
problems.

## Details

A diagnostic is named in `sampler_diagnostics_failed` if an R-hat value
exceeds 1.05 or if there is a divergent transition. It is named in
`sampler_diagnostics_warned` if an R-hat value lies between 1.01 and
1.05, a bulk or tail effective sample size is below 400, an iteration
reached the maximum tree depth, or the E-BFMI of a chain is below 0.2.
R-hat below 1.01 and effective sample sizes above 400 are the
recommendations of Vehtari et al. (2021); the other thresholds are
choices of this package.

With fewer than 400 draws in all, R-hat and the effective sample sizes
are too unreliable to act on. Their values are still returned, but they
are not compared with the thresholds, and `diagnostic_skipped` holds
`"rhat_check"`, `"ess_bulk_check"` and `"ess_tail_check"` to say so.
Divergent transitions, the maximum tree depth and E-BFMI are checked
however short the run.

R-hat and the effective sample sizes are undefined for a quantity that
takes the same value in every draw, such as the posterior mode of a site
whose mode never moves from one grid point, and are `NA` for it. If one
of the three diagnostics cannot be computed at all, or comes out
infinite for some quantity, all its values are `NA` and its name is in
`diagnostic_skipped`.
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
discusses what to do about each kind of warning.

## Functions

- `print(bef_diagnostic)`: Prints the worst value of each diagnostic.

- `summary(bef_diagnostic)`: Returns a list in which `rhat`, `ess_bulk`
  and `ess_tail` are reduced to their worst value, with the name of the
  quantity that has it.

- `format(bef_diagnostic)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows, as a character
  vector.

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
diag <- diagnose(raudenbush_fit)
diag
#> <bef_diagnostic>
#> Model family: RE
#> Rhat max: 1.007 (lp__)
#> ESS bulk min: 676.6 (lp__)
#> ESS tail min: 888.8 (lp__)
#> Divergences: 0
#> Max treedepth hits: 1
#> E-BFMI min: 0.2565
#> Runtime: 1.351 sec
#> Warnings: max_treedepth
#> Effective params: mean 10.29; sd 2.84

# The worst value of R-hat and the quantity that has it.
summary(diag)$rhat
#> $value
#> [1] 1.00688
#> 
#> $index
#> [1] 1
#> 
#> $variable
#> [1] "lp__"
#> 

# The three quantities with the smallest bulk effective sample size.
head(sort(diag$ess_bulk), 3)
#>     lp__   lambda log_g[8] 
#> 676.5790 808.8341 920.3033 
```
