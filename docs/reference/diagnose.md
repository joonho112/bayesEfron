# Extract sampler diagnostics from a bayesEfron fit

Return a structured `bef_diagnostic` object that bundles the
sampler-health and model-quality summaries already computed by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).
The function does not re-process draws or run new computations; it
surfaces information that the fit pipeline stored as attributes of
`fit$metadata`, behind the supported access path.

Use `diagnose()` together with `summary(fit)` to triage a fit:
[`summary()`](https://rdrr.io/r/base/summary.html) gives the user-facing
posterior summaries; `diagnose()` gives the structured sampler-health
view that pairs with `plot(fit, type = "diagnostic")`.

## Usage

``` r
diagnose(fit, ...)

# Default S3 method
diagnose(fit, ...)

# S3 method for class 'bef_fit'
diagnose(fit, ...)
```

## Arguments

- fit:

  A `bef_fit` object returned by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- ...:

  Reserved for future expansion; must be empty in v0.1.

## Value

A validated `bef_diagnostic` object with the schema tabulated in
**Details**.

## The `bef_diagnostic` schema

A `bef_diagnostic` object carries:

|  |  |  |
|----|----|----|
| Field | Type | Meaning |
| `rhat` | numeric vector | Per-parameter R-hat (Gelman–Rubin) statistic. Values close to 1 (typically `< 1.01`) indicate convergence. |
| `ess_bulk` | numeric vector | Per-parameter bulk effective sample size; large values support reliable posterior means. |
| `ess_tail` | numeric vector | Per-parameter tail ESS; large values support reliable tail-quantile / interval estimates. |
| `divergences` | numeric vector or NA | Number of divergent transitions per chain. |
| `max_treedepth` | numeric vector or NA | Number of max-treedepth saturations per chain. |
| `effective_params_summary` | named list | Posterior summary of effective parameters. |
| `model_family` | character | `"RE"` for v0.1. |
| `stan_file_sha256` | character | SHA-256 of the locked Stan source used for the fit. |
| `runtime_seconds` | numeric | Sampler wall-clock for the fit. |
| `diagnostic_skipped` | character vector | Diagnostics intentionally not computed (rare). |
| `sampler_diagnostics_failed` | character vector | Diagnostics requested but failed at extract time (rare). |

The returned object is validated internally before return.

## Reading the result

Call [`summary()`](https://rdrr.io/r/base/summary.html) on a
`bef_diagnostic` for the most-extreme value per field (max R-hat, min
ESS, total divergences, total max-treedepth saturations); call
[`print()`](https://rdrr.io/r/base/print.html) for a one-screen textual
view. The `plot(fit, type = "diagnostic")` view consumes the same
attribute payload through this generic.

## See also

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  for the upstream pipeline that populates these diagnostics.

- [summary.bef_diagnostic()](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
  for the extreme-value summary.

- [`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  for the graphical companion view.

- The methodological vignette M6 ("Verification and calibration") for
  what each diagnostic means in the context of the package's
  verification ledger.

## Examples

``` r
# Load the cached five-site smoke fit shipped with the package.
fit <- readRDS(system.file(
  "examples", "cached_fit_re_smoke.rds",
  package = "bayesEfron"
))

diag <- diagnose(fit)
print(diag)
#> <bef_diagnostic>
#> Model family: RE
#> Rhat max: 1.203
#> ESS bulk min: 12.53
#> ESS tail min: 105
#> Divergences: 0
#> Max treedepth hits: 0
#> Runtime: 0.2906 sec
#> Stan SHA-256: 57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e
#> Effective params: mean 0.3754; sd 0.3722
summary(diag)
#> $rhat
#> $rhat$value
#> [1] 1.203085
#> 
#> $rhat$index
#> [1] 1
#> 
#> 
#> $ess_bulk
#> $ess_bulk$value
#> [1] 12.53436
#> 
#> $ess_bulk$index
#> [1] 1
#> 
#> 
#> $ess_tail
#> $ess_tail$value
#> [1] 104.9594
#> 
#> $ess_tail$index
#> [1] 1
#> 
#> 
#> $divergences
#> [1] 0
#> 
#> $max_treedepth
#> [1] 0
#> 
#> $effective_params
#> $effective_params$mean
#> [1] 0.3753726
#> 
#> $effective_params$sd
#> [1] 0.3722022
#> 
#> $effective_params$q5
#> [1] 0.07486369
#> 
#> $effective_params$q50
#> [1] 0.3212756
#> 
#> $effective_params$q95
#> [1] 0.7056425
#> 
#> 
#> $model_family
#> [1] "RE"
#> 
#> $stan_file_sha256
#> [1] "57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e"
#> 
#> $runtime_seconds
#> [1] 0.2905509
#> 
#> $diagnostic_skipped
#> character(0)
#> 
#> $sampler_diagnostics_failed
#> character(0)
#> 
```
