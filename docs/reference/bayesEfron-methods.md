# S3 methods for bayesEfron objects

The package exports a coordinated set of S3 methods for the `bef_fit`,
`bef_fit_re`, `bef_data`, and `bef_diagnostic` classes. Together they
cover printing, summarising, extracting credible intervals, pulling
point estimates, converting to a data frame, counting observations,
computing the marginal log-likelihood, and coercing to the
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
format.

This help page documents the family-agnostic surface, which works on any
fitted `bef_fit` object regardless of the model family. The RE-specific
child class adds [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and a
refined [`summary()`](https://rdrr.io/r/base/summary.html). Plotting is
documented separately at
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md);
the diagnostic producer is at
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md).

## Usage

``` r
# S3 method for class 'bef_fit'
format(x, ..., use_cli = NULL)

# S3 method for class 'summary.bef_fit'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_data'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_diagnostic'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_fit'
print(x, ...)

# S3 method for class 'bef_fit'
summary(object, level = 0.9, ...)

# S3 method for class 'summary.bef_fit'
print(x, ...)

# S3 method for class 'bef_fit_re'
summary(object, level = 0.9, ...)

# S3 method for class 'bef_fit_re'
coef(object, type = c("mean", "map"), ...)

# S3 method for class 'bef_fit_re'
vcov(object, ...)

# S3 method for class 'bef_fit_re'
confint(object, parm = NULL, level = 0.9, type = c("theta", "g"), ...)

# S3 method for class 'bef_fit_re'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'bef_data'
print(x, ...)

# S3 method for class 'bef_data'
summary(object, ...)

# S3 method for class 'bef_diagnostic'
print(x, ...)

# S3 method for class 'bef_diagnostic'
summary(object, ...)

# S3 method for class 'bef_fit'
nobs(object, ...)

# S3 method for class 'bef_fit'
logLik(object, ...)

# S3 method for class 'bef_fit'
as_draws(x, ...)
```

## Arguments

- x, object:

  A bayesEfron S3 object (`bef_fit`, `bef_fit_re`, `bef_data`, or
  `bef_diagnostic`, depending on the method).

- ...:

  Additional arguments passed to methods. Most methods ignore `...`;
  [`format()`](https://rdrr.io/r/base/format.html)-based printers accept
  it for compatibility with the generic.

- use_cli:

  `NULL`, `TRUE`, or `FALSE`; controls optional `cli`-styled output for
  [`format()`](https://rdrr.io/r/base/format.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods. See the "Print
  and format backends" subsection in **Details**.

- level:

  Numeric credible level in \\(0, 1)\\. Defaults to `0.9` (90 percent
  credible interval).

- type:

  Character option for methods with multiple estimands:
  [`coef()`](https://rdrr.io/r/stats/coef.html) accepts `"mean"`
  (posterior mean, default) or `"map"` (posterior mode);
  [`confint()`](https://rdrr.io/r/stats/confint.html) accepts `"theta"`
  (per-site latent-effect intervals, default) or `"g"`
  (mixing-distribution functionals).

- parm:

  Optional parameter or site subset. `NULL` (default) returns all sites
  or all `g`-functionals. Numeric values are 1-based site indices;
  character values are matched against the site labels.

- row.names:

  Optional character vector of row names for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), with
  one value per site.

- optional:

  Included for
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
  compatibility; not used.

## Value

The return type depends on the method; see the index table in
**Details**.

## Method index

|  |  |  |
|----|----|----|
| Method | Class dispatched on | Returns |
| [`print()`](https://rdrr.io/r/base/print.html) | `bef_fit`, `bef_data`, `bef_diagnostic`, `summary.bef_fit` | invisibly the input |
| [`summary()`](https://rdrr.io/r/base/summary.html) | `bef_fit`, `bef_fit_re`, `bef_data`, `bef_diagnostic` | a list (class `summary.bef_fit*` for fits) |
| [`format()`](https://rdrr.io/r/base/format.html) | `bef_fit`, `summary.bef_fit`, `bef_data`, `bef_diagnostic` | character vector |
| [`coef()`](https://rdrr.io/r/stats/coef.html) | `bef_fit_re` | named numeric vector of site point estimates |
| [`vcov()`](https://rdrr.io/r/stats/vcov.html) | `bef_fit_re` | diagonal matrix of site posterior variances |
| [`confint()`](https://rdrr.io/r/stats/confint.html) | `bef_fit_re` | data frame of credible intervals |
| [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) | `bef_fit_re` | the per-site `theta_summary` table |
| [`nobs()`](https://rdrr.io/r/stats/nobs.html) | `bef_fit` | integer site count `K` |
| [`logLik()`](https://rdrr.io/r/stats/logLik.html) | `bef_fit` | `logLik` object with `df` = effective parameters |
| [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html) | `bef_fit` | `draws_array` |

## Metadata access

`fit$metadata` is a closed list of 13 named fields; four additional
payloads are stored as attributes (`diagnostics`, `diagnostic_skipped`,
`sampler_diagnostics_failed`, `sd_g_summary`). The methods on this page
are the recommended access path:
[`summary()`](https://rdrr.io/r/base/summary.html) surfaces the
user-relevant fields;
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
returns a structured `bef_diagnostic` object covering the diagnostic
attributes; direct [`attr()`](https://rdrr.io/r/base/attr.html) access
is supported but unnecessary.

## Sites and labels

Site labels propagated from the input (named list or `escalc` row names)
appear in [`coef()`](https://rdrr.io/r/stats/coef.html) names,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) dimnames, the `site`
column of [`confint()`](https://rdrr.io/r/stats/confint.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
[`summary()`](https://rdrr.io/r/base/summary.html) print output. When
labels are absent the numeric site index is used.

## Print and format backends

[`print()`](https://rdrr.io/r/base/print.html) and
[`format()`](https://rdrr.io/r/base/format.html) for `bef_fit`,
`summary.bef_fit`, `bef_data`, and `bef_diagnostic` prefer `cli`-styled
output when the `cli` package is installed and stdout supports it, and
fall back to a plain base-R representation otherwise. The `use_cli`
argument lets callers force one branch:

- `use_cli = NULL` (default) — auto-detect.

- `use_cli = TRUE` — force `cli` styling (errors if `cli` is not
  installed).

- `use_cli = FALSE` — force plain base output.

Setting the environment variable `BAYESEFRON_NO_CLI=1` is equivalent to
`use_cli = FALSE` for every call in the session and is the recommended
way to suppress styling in CI logs and redirected stdout.

## Site count requirement

The `bef_data` validator (used by
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
and
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md))
requires at least five sites. The constraint reflects the methodological
requirement of the log-spline deconvolution prior, not a software
limitation.

## See also

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  for producing the fitted object.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  for the structured diagnostic producer.

- [`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  for visualization.

## Examples

``` r
# Load the cached five-site smoke fit shipped with the package.
fit <- readRDS(system.file(
  "examples", "cached_fit_re_smoke.rds",
  package = "bayesEfron"
))

summary(fit)
#> <summary.bef_fit>
#> 
#> Prior g:
#>   mean: 0.174
#>   var:  0.2094
#>   sd:   0.4549
#> 
#> Diagnostics:
#>   Rhat:              1.203
#>   ESS bulk:          12.53
#>   ESS tail:          105
#>   Divergences:       0
#>   Max treedepth:     0
#>   Effective params:  mean 0.3754, sd 0.3722
#>   Log marginal lik.: mean -2.433, sd 0.443
#>   Runtime:           0.2906 sec
#>   Stan SHA-256:      57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e
#> 
#> Theta summary:
#>   site      mean        sd     lower     upper       map
#>   1      -0.1948    0.1713    -0.456     0.036   -0.1956
#>   2       0.0401    0.1471   -0.1608       0.2   0.03862
#>   3       0.1808    0.2128   -0.1952    0.4968    0.1787
#>   4       0.3652    0.1869   0.06716    0.6608     0.363
#>   5       0.5658    0.2197    0.2984    0.9216    0.5847
print(fit)
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 5
#> Grid method: paper_realdata
#> Runtime: 0.2906 sec
#> Stan SHA-256: 57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e
#> Diagnostics: Rhat 1.203; ESS bulk 12.53; ESS tail 105; divergences 0; max treedepth 0
#> Use summary() for posterior summaries.
coef(fit)                       # posterior mean per site
#>           1           2           3           4           5 
#> -0.19484562  0.04010493  0.18079325  0.36515025  0.56577011 
coef(fit, type = "map")         # posterior MAP per site
#>         1         2         3         4         5 
#> -0.195568  0.038624  0.178680  0.363016  0.584744 
confint(fit)                    # 90 percent credible intervals on theta
#>   site    lower   upper       point
#> 1    1 -0.45600 0.03600 -0.19484562
#> 2    2 -0.16080 0.20000  0.04010493
#> 3    3 -0.19524 0.49684  0.18079325
#> 4    4  0.06716 0.66084  0.36515025
#> 5    5  0.29840 0.92160  0.56577011
confint(fit, level = 0.95, type = "g")
#>     site      lower     upper     point
#> 1 mean_g 0.04313042 0.2983847 0.1739745
#> 2  var_g 0.10362327 0.2567256 0.2094317
#> 3   sd_g 0.31627393 0.5066809 0.4548573
nobs(fit)                       # site count K
#> [1] 5
logLik(fit)
#> 'log Lik.' -2.432978 (df=0.3753726)
as.data.frame(fit)              # per-site summary as a data frame
#>   site        mean        sd hpdi_lower hpdi_upper       map
#> 1    1 -0.19484562 0.1713360   -0.45600    0.03600 -0.195568
#> 2    2  0.04010493 0.1471210   -0.16080    0.20000  0.038624
#> 3    3  0.18079325 0.2128270   -0.19524    0.49684  0.178680
#> 4    4  0.36515025 0.1868645    0.06716    0.66084  0.363016
#> 5    5  0.56577011 0.2196747    0.29840    0.92160  0.584744
posterior::as_draws(fit)        # draws_array for downstream tools
#> # A draws_array: 100 iterations, 1 chains, and 183 variables
#> , , variable = lp__
#> 
#>          chain
#> iteration     1
#>         1  0.14
#>         2 -1.41
#>         3 -3.43
#>         4 -2.63
#>         5 -5.09
#> 
#> , , variable = alpha[1]
#> 
#>          chain
#> iteration     1
#>         1 -0.32
#>         2  0.61
#>         3  0.94
#>         4  0.95
#>         5  0.66
#> 
#> , , variable = alpha[2]
#> 
#>          chain
#> iteration     1
#>         1 -0.29
#>         2  0.81
#>         3  1.28
#>         4  1.10
#>         5  0.41
#> 
#> , , variable = alpha[3]
#> 
#>          chain
#> iteration      1
#>         1  0.148
#>         2 -0.081
#>         3 -0.023
#>         4 -0.139
#>         5  1.010
#> 
#> # ... with 95 more iterations, and 179 more variables
```
