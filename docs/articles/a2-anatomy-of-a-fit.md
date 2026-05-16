# A2 · Anatomy of a bayesEfron fit

## Scope

The **A1 · Getting started** vignette ran a five-site fit, called
[`summary()`](https://rdrr.io/r/base/summary.html), and read off the two
inferential targets — the mixing distribution $`g(\theta)`$ and the
per-site posteriors $`\theta_i \mid
\hat\theta_i, \sigma_i`$. This vignette opens up the same object. Every
visible field of the returned `bef_fit_re`, every attribute payload, and
every S3 method is named, shown, and explained.

All examples below use the cached five-site smoke fit shipped with the
package:

``` r

fit <- readRDS(system.file(
  "examples", "cached_fit_re_smoke.rds",
  package = "bayesEfron"
))
```

No code chunk in this vignette compiles or samples from Stan; the
fixture stands in for a live
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
call.

## Object structure

The returned object carries two S3 classes — a family-agnostic parent
and an RE-specific child — and three top-level slots.

``` r

class(fit)
#> [1] "bef_fit_re" "bef_fit"
```

``` r

str(fit, max.level = 1)
#> List of 3
#>  $ draws    : 'draws_array' num [1:100, 1, 1:183] 0.141 -1.412 -3.427 -2.627 -5.091 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ metadata :List of 13
#>   ..- attr(*, "sd_g_summary")=List of 5
#>   ..- attr(*, "diagnostics")=List of 5
#>   ..- attr(*, "diagnostic_skipped")= chr(0) 
#>   ..- attr(*, "sampler_diagnostics_failed")= chr(0) 
#>  $ posterior:List of 9
#>  - attr(*, "class")= chr [1:2] "bef_fit_re" "bef_fit"
```

The three top-level fields are coordinated. `draws` is the raw
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
of MCMC draws across all monitored parameters and generated quantities.
`metadata` is a closed list of exactly 13 named fields, with four
additional payloads stored as attributes. `posterior` is a tidy,
by-parameter representation of the draws that the S3 methods consume; it
does not replace `draws`, but flattens what those methods need into a
faster-to-index form.

## The 13 visible metadata fields

The metadata contract is closed:
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
always populates exactly the same 13 names in the same order. Treat any
departure as a bug.

``` r

names(fit$metadata)
#>  [1] "model_family"                    "grid_method"                    
#>  [3] "seed"                            "cmdstan_version"                
#>  [5] "stan_file_sha256"                "data_list"                      
#>  [7] "runtime_seconds"                 "mean_g_summary"                 
#>  [9] "var_g_summary"                   "theta_summary"                  
#> [11] "theta_rep_draws"                 "effective_params_summary"       
#> [13] "log_marginal_likelihood_summary"
```

The 13 fields divide naturally into five groupings.

| Field | Purpose |
|:---|:---|
| `model_family` | Model family identifier; `"RE"` only in v0.1. |
| `grid_method` | Which of the four grid recipes was used. |
| `seed`, `cmdstan_version`, `stan_file_sha256` | Reproducibility triple: effective seed, CmdStan version, SHA-256 of the locked Stan source. |
| `data_list` | The seven-field Stan data block sent to CmdStan. |
| `runtime_seconds` | Sampler wall-clock in seconds. |
| `mean_g_summary`, `var_g_summary` | Posterior summaries of two functionals of the mixing distribution $`g`$. |
| `theta_summary` | Per-site posterior table: one row per site, columns include `mean`, `sd`, `hpdi_lower`, `hpdi_upper`, `map`. |
| `theta_rep_draws` | Replicated effects ($`J \times K`$) used for the posterior-predictive `confint(type = "theta")` path. |
| `effective_params_summary`, `log_marginal_likelihood_summary` | Model-quality summaries: effective parameter count and log marginal likelihood. |

The 13 visible metadata fields of `bef_fit_re`, grouped by purpose.

The reproducibility triple is the part that matters for archival. These
fields, together with `data_list`, sampler settings, package version,
and a matching CmdStan/toolchain environment, make the fit auditable and
practically replayable. They should not be described as a cross-platform
guarantee of byte-identical draws.

`data_list` is a passthrough record of the Stan data block, exposed so
that downstream tools can recover the exact inputs CmdStan saw without
re-running the grid recipe:

``` r

names(fit$metadata$data_list)
#> [1] "K"         "theta_hat" "sigma"     "L"         "grid"      "M"        
#> [7] "B"
```

`theta_summary` is the per-site table that
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
`confint(type = "theta")` ultimately read from; `theta_rep_draws` is the
matrix of posterior-predictive replicates that powers
`confint(type = "theta")`. Both are documented as part of the metadata
contract precisely so users can reach for them directly when an
unsupported aggregation is needed.

## The four attribute payloads

Four further objects live as attributes of `fit$metadata` rather than as
named slots: one secondary functional of $`g`$, one sampler-diagnostic
record, and two diagnostic bookkeeping flags.

``` r

setdiff(names(attributes(fit$metadata)), "names")
#> [1] "sd_g_summary"               "diagnostics"               
#> [3] "diagnostic_skipped"         "sampler_diagnostics_failed"
```

The package treats these as durable but not part of the print surface.
They are read through dedicated S3 methods rather than direct
[`attr()`](https://rdrr.io/r/base/attr.html) calls:

- `sd_g_summary` — posterior summary of $`\mathrm{sd}(g)`$; surfaced by
  [`summary()`](https://rdrr.io/r/base/summary.html) (as
  `prior_summary$sd`) and by `confint(type = "g")`.
- `diagnostics` — the HMC diagnostic record: `rhat`, `ess_bulk`,
  `ess_tail`, `divergences`, `max_treedepth`. Surfaced by
  [`summary()`](https://rdrr.io/r/base/summary.html) and by
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md).
- `diagnostic_skipped` — character vector of diagnostic names that the
  post-processor refused to compute (empty on a clean fit).
- `sampler_diagnostics_failed` — character vector of HMC sampler
  diagnostics that the underlying CmdStanMCMC call could not return
  (empty on a clean fit).

For the cached smoke fit both bookkeeping flags are empty character
vectors, which is the expected shape on any clean run:

``` r

attr(fit$metadata, "diagnostic_skipped")
#> character(0)
attr(fit$metadata, "sampler_diagnostics_failed")
#> character(0)
```

A non-empty value in either flag is the post-processor’s way of
recording that one or more diagnostic computations could not complete;
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
propagate the names without guessing at substitute values.

## The S3 method surface

The methods table from `?bayesEfron-methods` enumerates the coordinated
method set. The subsections below show one concrete call per method
against the cached fit.

### `summary(fit)` and `print(fit)`

`summary(fit)` is the user-facing entry point: it returns a three-block
`summary.bef_fit_re` object covering the prior $`g`$ summary, the HMC
diagnostic record, and the per-site `theta_summary` table. Its output
was the subject of the **A1** walkthrough and is not repeated here.

`print(fit)` returns the compact one-screen overview that fits in a
console snippet. It does not enumerate sites; for that, call
[`summary()`](https://rdrr.io/r/base/summary.html) or
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

``` r

print(fit)
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 5
#> Grid method: paper_realdata
#> Runtime: 0.2906 sec
#> Stan SHA-256: 57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e
#> Diagnostics: Rhat 1.203; ESS bulk 12.53; ESS tail 105; divergences 0; max treedepth 0
#> Use summary() for posterior summaries.
```

### `coef(fit)`: point estimates

[`coef()`](https://rdrr.io/r/stats/coef.html) dispatches on the
RE-specific child class and returns a named numeric vector of per-site
point estimates. The default, `type = "mean"`, returns posterior means.

``` r

coef(fit)
#>           1           2           3           4           5 
#> -0.19484562  0.04010493  0.18079325  0.36515025  0.56577011
```

`type = "map"` returns the posterior mode (maximum a posteriori) column
from the same table.

``` r

coef(fit, type = "map")
#>         1         2         3         4         5 
#> -0.195568  0.038624  0.178680  0.363016  0.584744
```

The two columns rarely coincide exactly, but for a unimodal per-site
posterior on a continuous support they typically agree to within a
fraction of $`\sigma_i`$.

### `vcov(fit)`: diagonal posterior covariance

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns a $`K \times K`$
diagonal matrix whose entries are the squared posterior standard
deviations of the per-site $`\theta_i`$. The matrix is diagonal by
construction: the package does not model posterior dependence across
sites at the $`\theta_i`$ level.

``` r

vcov(fit)
#>            1         2          3          4          5
#> 1 0.02935603 0.0000000 0.00000000 0.00000000 0.00000000
#> 2 0.00000000 0.0216446 0.00000000 0.00000000 0.00000000
#> 3 0.00000000 0.0000000 0.04529531 0.00000000 0.00000000
#> 4 0.00000000 0.0000000 0.00000000 0.03491835 0.00000000
#> 5 0.00000000 0.0000000 0.00000000 0.00000000 0.04825697
```

The dimnames are the site labels (or the numeric site indices when
labels are absent, as in this smoke fixture).

### `confint(fit)` and `confint(fit, type = "g")`

`confint(fit)` returns a data frame of credible intervals on the
per-site latent effects $`\theta_i`$. The default credible level is
$`90\%`$ (`level = 0.9`); the intervals are read directly from the
`theta_rep_draws` posterior-predictive matrix.

``` r

confint(fit)
#>   site    lower   upper       point
#> 1    1 -0.45600 0.03600 -0.19484562
#> 2    2 -0.16080 0.20000  0.04010493
#> 3    3 -0.19524 0.49684  0.18079325
#> 4    4  0.06716 0.66084  0.36515025
#> 5    5  0.29840 0.92160  0.56577011
```

Switching `type = "g"` returns intervals on the three functionals of the
mixing distribution stored under `posterior`: the mean of $`g`$, the
variance of $`g`$, and the standard deviation of $`g`$. The `level`
argument controls the credible width independently of the type.

``` r

confint(fit, level = 0.95, type = "g")
#>     site      lower     upper     point
#> 1 mean_g 0.04313042 0.2983847 0.1739745
#> 2  var_g 0.10362327 0.2567256 0.2094317
#> 3   sd_g 0.31627393 0.5066809 0.4548573
```

The `parm` argument (not shown here) accepts numeric site indices,
character site labels, or character functional names (`"mean_g"`,
`"var_g"`, `"sd_g"`) for selecting a subset.

### `nobs(fit)`: site count

[`nobs()`](https://rdrr.io/r/stats/nobs.html) is dispatched on the
parent class and returns the integer site count $`K`$ pulled from
`fit$metadata$data_list$K`.

``` r

nobs(fit)
#> [1] 5
```

### `logLik(fit)`: log marginal likelihood

[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns a base-R
`logLik` object whose value is the posterior mean of the log marginal
likelihood. The `df` attribute holds the posterior mean of the effective
parameter count, and `nobs` is the site count. The object is exposed for
inspection and for user-defined comparisons; formal AIC/BIC or LOO-CV
guidance is not part of v0.1.0.

``` r

logLik(fit)
#> 'log Lik.' -2.432978 (df=0.3753726)
```

### `as.data.frame(fit)`: per-site table

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `theta_summary` table directly, with one row per site and the same
`site`, `mean`, `sd`, `hpdi_lower`, `hpdi_upper`, and `map` columns that
[`summary()`](https://rdrr.io/r/base/summary.html) surfaces.

``` r

as.data.frame(fit)
#>   site        mean        sd hpdi_lower hpdi_upper       map
#> 1    1 -0.19484562 0.1713360   -0.45600    0.03600 -0.195568
#> 2    2  0.04010493 0.1471210   -0.16080    0.20000  0.038624
#> 3    3  0.18079325 0.2128270   -0.19524    0.49684  0.178680
#> 4    4  0.36515025 0.1868645    0.06716    0.66084  0.363016
#> 5    5  0.56577011 0.2196747    0.29840    0.92160  0.584744
```

This is the method to reach for when feeding a `bayesEfron` fit into a
downstream pipeline that expects a tabular per-site object; it preserves
the original column types and is stable across releases.

### `posterior::as_draws(fit)`: draws array

For workflows that want to drive the draws directly through the
`posterior` or `bayesplot` ecosystems, the
[`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html)
generic returns the underlying `draws_array`.

``` r

da <- posterior::as_draws(fit)
class(da)
#> [1] "draws_array" "draws"       "array"
dim(da)
#> [1] 100   1 183
```

The array dimensions are `iterations × chains × variables`. For the
cached smoke fit this is $`100 \times 1 \times 183`$ — 100 sampling
iterations from a single chain, with 183 monitored parameters and
generated quantities.

## Diagnostics object

[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
is the structured producer for the `bef_diagnostic` class. It pulls the
`diagnostics` attribute together with the durable model-quality
summaries and emits a single object that can be passed around, printed,
and inspected without re-reaching into the metadata attributes.

``` r

diag <- diagnose(fit)
class(diag)
#> [1] "bef_diagnostic"
diag
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
```

The printed view shows the worst-case $`\widehat{R}`$, the smallest bulk
and tail ESS, the total divergence count, and the number of
maximum-treedepth hits, alongside the runtime, the Stan SHA-256, and the
effective parameter summary. The full diagnostic checklist — including
release thresholds and remediation patterns — is the subject of the **A5
· Diagnostics in practice** vignette.

## What’s next

The remaining release vignettes carry the workflow forward along four
orthogonal axes.

- **A3 · From
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  to bayesEfron** — the end-to-end path that begins with raw treatment /
  control group summaries, computes effect sizes with
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html),
  and lands at a fitted object whose anatomy is the one shown here.
- **A4 · Choosing a grid recipe** — when to use `"paper_realdata"`, when
  an oracle-aware `"paper_simulation"` or `"paper_sensitivity"` run is
  appropriate, and what the experimental `"kl_target_experimental"`
  recipe actually targets.
- **A5 · Diagnostics in practice** — the release-quality checklist for
  the `diagnostics` attribute and the `bef_diagnostic` object, including
  the conventional $`\widehat{R}`$ and ESS thresholds and how to act on
  divergence or treedepth flags.
- **A6 · Plotting `bayesEfron` fits** — the
  [`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  surface, including the caterpillar plot used in **A1**, the
  prior-$`g`$ density, and the posterior-predictive plots that read off
  `theta_rep_draws`.
