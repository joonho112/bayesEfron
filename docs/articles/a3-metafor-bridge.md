# A3 · From metafor::escalc() to bayesEfron

## Scope

The **A1 · Getting started** vignette began at the effect-size scale: a
length-five `theta_hat` vector and a matching `sigma` vector dropped
straight into
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).
Many applied meta-analyses do not begin there. The raw record is per-arm
summaries — group means, standard deviations, and sample sizes — and the
analyst’s first task is to convert those summaries into a per-site
estimate and its sampling standard error. This vignette walks the full
path: from the raw aggregate data, through
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html),
through the
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
bridge, into
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
and out to the standard outputs.

## The starting point

Consider a randomized education trial run across $`K = 10`$ schools.
Each school contributes a treatment arm and a control arm, with each arm
reporting an outcome mean, the within-arm standard deviation, and the
per-arm student count. The outcome is a standardized end-of-year test
score; the program assignment was randomized within each school. The
substantive research question has two parts. First, what does the
*population* distribution of school-level program effects look like —
concentrated around a common effect, or genuinely heterogeneous? Second,
what is each *individual* school’s program effect after partial pooling
absorbs the noise in its small-sample estimate?

The raw record is the following data frame. Column names follow the
convention used by
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
for two-sample summaries: `m1i`, `sd1i`, `n1i` for the treatment arm and
`m2i`, `sd2i`, `n2i` for the control arm.

``` r

trial <- data.frame(
  site = paste0("site_", 1:10),
  m1i  = c(0.42, 0.18, 0.05, 0.31, 0.54,
           0.12, 0.39, 0.27, 0.48, 0.16),
  sd1i = c(1.05, 1.12, 0.98, 1.08, 1.15,
           1.03, 1.10, 0.96, 1.20, 1.07),
  n1i  = c(42L, 39L, 44L, 41L, 40L,
           43L, 37L, 45L, 38L, 46L),
  m2i  = c(0.05, 0.02, -0.08, 0.04, 0.10,
           -0.03, 0.11, 0.00, 0.08, -0.02),
  sd2i = c(1.01, 1.09, 1.04, 1.02, 1.12,
           1.00, 1.06, 0.99, 1.16, 1.05),
  n2i  = c(40L, 38L, 42L, 39L, 41L,
           44L, 36L, 43L, 39L, 45L)
)
trial
#>       site  m1i sd1i n1i   m2i sd2i n2i
#> 1   site_1 0.42 1.05  42  0.05 1.01  40
#> 2   site_2 0.18 1.12  39  0.02 1.09  38
#> 3   site_3 0.05 0.98  44 -0.08 1.04  42
#> 4   site_4 0.31 1.08  41  0.04 1.02  39
#> 5   site_5 0.54 1.15  40  0.10 1.12  41
#> 6   site_6 0.12 1.03  43 -0.03 1.00  44
#> 7   site_7 0.39 1.10  37  0.11 1.06  36
#> 8   site_8 0.27 0.96  45  0.00 0.99  43
#> 9   site_9 0.48 1.20  38  0.08 1.16  39
#> 10 site_10 0.16 1.07  46 -0.02 1.05  45
```

Treatment-arm means range from $`0.05`$ to $`0.54`$, control-arm means
from $`-0.08`$ to $`0.11`$, and per-arm sample sizes from $`36`$ to
$`46`$ students. Three features of the record matter for the bridge. The
two arms at each site are independent random samples, so the sampling
variance of the per-site mean difference is the sum of the two per-arm
sampling variances;
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
computes this for the analyst. The outcome scale is the standardized
test score, so the mean-difference effect size is interpretable on the
same scale as a Cohen’s-$`d`$-style standardized mean difference without
an additional pooling step. The per-arm sample sizes are modest enough
that within-site sampling noise is large relative to plausible true
between-site variation — the regime in which a hierarchical estimator
earns its keep.

## Effect sizes via `metafor::escalc()`

The `metafor` package is the canonical R interface for the family of
per-study effect-size and sampling-variance computations that upstream a
meta-analysis. The `escalc()` constructor takes the group-arm summaries
and a `measure` argument that selects the estimator. For continuous
outcomes on a common scale `measure = "MD"` (the raw mean difference) is
the natural choice; for studies that report on different outcome scales
`measure = "SMD"` (the standardized mean difference, Hedges’s $`g`$)
would normalize them to a common standard-deviation unit. The trial
above uses a single standardized test, so the raw mean difference is
what we want.

``` r

esc <- metafor::escalc(
  measure = "MD",
  m1i  = m1i,  sd1i = sd1i, n1i = n1i,
  m2i  = m2i,  sd2i = sd2i, n2i = n2i,
  data = trial,
  slab = site
)
class(esc)
#> [1] "escalc"     "data.frame"
```

The returned object has class `c("escalc", "data.frame")` and carries
two new columns alongside the original ones.

``` r

knitr::kable(
  as.data.frame(esc)[, c("site", "yi", "vi")],
  digits  = 4,
  caption = paste(
    "Per-site effect sizes from `metafor::escalc(measure = \"MD\")`.",
    "Column `yi` is the treatment-minus-control mean difference;",
    "column `vi` is its sampling variance, computed as the sum of",
    "the two per-arm variances of the mean."
  )
)
```

| site    |   yi |     vi |
|:--------|-----:|-------:|
| site_1  | 0.37 | 0.0518 |
| site_2  | 0.16 | 0.0634 |
| site_3  | 0.13 | 0.0476 |
| site_4  | 0.27 | 0.0551 |
| site_5  | 0.44 | 0.0637 |
| site_6  | 0.15 | 0.0474 |
| site_7  | 0.28 | 0.0639 |
| site_8  | 0.27 | 0.0433 |
| site_9  | 0.40 | 0.0724 |
| site_10 | 0.18 | 0.0494 |

Per-site effect sizes from `metafor::escalc(measure = "MD")`. Column
`yi` is the treatment-minus-control mean difference; column `vi` is its
sampling variance, computed as the sum of the two per-arm variances of
the mean.

Read the two new columns directly. `yi` is the per-site estimate
$`\widehat\theta_i = \bar y_{i,T} - \bar y_{i,C}`$: the
treatment-minus-control mean difference for site $`i`$. `vi` is the
corresponding sampling variance
$`\widehat{\mathrm{Var}}(\widehat\theta_i) =
s_{i,T}^2/n_{i,T} + s_{i,C}^2/n_{i,C}`$: the textbook estimate of the
sampling variance of an independent-samples mean difference. The
sampling *standard error* — the quantity `bayesEfron` needs — is
$`\sigma_i = \sqrt{v_i}`$. The `slab` argument records the site labels
as an attribute of `yi`; they will travel into `escalc`’s print method
and into downstream `metafor` models, but `escalc` itself does not
promote them to the data frame’s row names.

The per-site estimates from this trial run from $`0.13`$ to $`0.44`$ and
the standard errors from roughly $`0.21`$ to $`0.27`$, all on the
standardized-score scale. The classical fixed-effect summary would be
the inverse-variance-weighted mean; the random-effects summary would be
a method-of-moments or REML estimate of the between-site variance. Both
throw away the per-site point estimates after a single pass.
`bayesEfron`’s contribution is to keep them: it estimates a continuous
mixing distribution $`g(\theta)`$ over the latent site effects and
returns a posterior for each $`\theta_i`$ separately, in addition to
summaries of $`g`$ itself.

## The `as_bef_data()` bridge

The
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
generic adapts the `escalc` object into the `bef_data` class that
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
consumes. The method extracts `theta_hat` from `yi`, computes `sigma` as
`sqrt(vi)` after a strict-positivity check on `vi`, and records a
provenance tag.

``` r

bef_dat <- as_bef_data(esc)
bef_dat
#> <bef_data>
#> Sites: 10
#> Source: metafor::escalc
#> theta_hat: min 0.13; median 0.27; max 0.44
#> sigma: min 0.208; median 0.2311; max 0.2691
```

The compact print shows the four pieces that survived the conversion:
the site count, the source tag, and one-line min / median / max
summaries of `theta_hat` and `sigma`. The full representation is a
four-element list:

``` r

names(bef_dat)
#> [1] "theta_hat" "sigma"     "names"     "source"
bef_dat$source
#> [1] "metafor::escalc"
```

The `source = "metafor::escalc"` element is the provenance tag the
`escalc` method writes. Downstream consumers can branch on it when the
bridge matters — for example, to record the input shape in a methods
appendix without re-deriving it from object class. `bef_dat$names` is
`NULL` here because the `escalc` object’s row names are the
integer-sequence default; the `slab` attribute that travels on `yi` is
the `metafor`-internal label, not a row name, and is dropped by the
bridge. Where the analyst wants the site labels surfaced in the fit’s
per-site table, the simplest fix is to set `row.names(esc) <- esc$site`
before the
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
call.

## Cached fixture

Running the full fit on the ten-site escalc data above requires a
working CmdStan installation and a non-trivial sampler run. So that this
vignette renders without a CmdStan dependency, the posterior-bearing
outputs below are read from the cached five-site smoke fixture shipped
at `inst/examples/cached_fit_re_smoke.rds`. The fixture was produced
with $`L = 51`$, $`M = 3`$, one chain, 150 warmup iterations, and 100
sampling iterations on the five-site toy data of **A1**. It is
deliberately small and is reserved for compile-and-run verification; the
diagnostic profile reported below reflects that budget, not a production
fit.

The displayed outputs —
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
the caterpillar plot — therefore correspond to the five-site smoke
configuration, not to the ten-site escalc data computed in the previous
section. The intent is to show the *shape* of each return value with a
renderable fixture; the **Live alternative** section below shows the
call you would make against the ten-site data when CmdStan is available
locally.

``` r

fit <- readRDS(system.file(
  "examples", "cached_fit_re_smoke.rds",
  package = "bayesEfron"
))
class(fit)
#> [1] "bef_fit_re" "bef_fit"
```

## Outputs

`summary(fit)` returns a three-block object with the posterior summary
of the mixing distribution $`g`$, the HMC diagnostic record, and the
per-site `theta_summary` table.

``` r

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
```

`confint(fit)` returns the per-site credible intervals on the latent
effects $`\theta_i`$. The default credible level is $`0.90`$ and the
default `type` is `"theta"`; the intervals are read from the
posterior-predictive `theta_rep_draws` matrix stored in metadata. The
same generic with `type = "g"` switches to credible intervals on the
functionals of the mixing distribution.

``` r

confint(fit, level = 0.9, type = "theta")
#>   site    lower   upper       point
#> 1    1 -0.45600 0.03600 -0.19484562
#> 2    2 -0.16080 0.20000  0.04010493
#> 3    3 -0.19524 0.49684  0.18079325
#> 4    4  0.06716 0.66084  0.36515025
#> 5    5  0.29840 0.92160  0.56577011
```

`diagnose(fit)` produces a `bef_diagnostic` object that bundles the HMC
sampler health record together with the durable model-quality summaries
(effective parameter count, log marginal likelihood). The print view
surfaces the worst-case $`\widehat{R}`$, the smallest bulk and tail ESS,
the divergence and maximum-treedepth counts, and the Stan SHA-256.
**A5** documents the release-quality checklist this object is built to
support.

``` r

diagnose(fit)
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

`plot(fit, type = "caterpillar")` is the visual summary of the per-site
posterior table. Each site contributes one point estimate and one
credible interval, ordered by the site index.

``` r

plot(fit, type = "caterpillar")
```

![Per-site posterior latent effects with 90 percent credible intervals
from the cached five-site smoke fit. The fixture is a smoke-scale
configuration (one chain, 100 sampling iterations) reserved for
compile-and-run verification; production fits with the package defaults
provide a larger post-warmup draw set and more reliable
diagnostics.](a3-metafor-bridge_files/figure-html/caterpillar-1.png)

Per-site posterior latent effects with 90 percent credible intervals
from the cached five-site smoke fit. The fixture is a smoke-scale
configuration (one chain, 100 sampling iterations) reserved for
compile-and-run verification; production fits with the package defaults
provide a larger post-warmup draw set and more reliable diagnostics.

## Live alternative

The chunk below is the call that produces a fitted object for the
ten-site escalc data computed earlier in this vignette. It runs at the
same smoke-scale configuration as the cached fixture and requires a
working CmdStan installation. Set the environment variable
`BAYESEFRON_RUN_LIVE=1` before rendering to execute it inline; otherwise
it is skipped.

``` r

fit_live <- bayes_efron_fit(
  theta_hat     = bef_dat$theta_hat,
  sigma         = bef_dat$sigma,
  L             = 51L,
  M             = 3L,
  chains        = 1L,
  iter_warmup   = 150L,
  iter_sampling = 100L,
  seed          = 1234L
)
summary(fit_live)
```

The same four S3 methods shown above —
[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) — apply
unchanged to `fit_live`. The canonical executable form of this
end-to-end path is installed at
`system.file("examples", "example-bayes-efron-fit.R", package = "bayesEfron")`
and is the script the vignette’s outputs would match in a live render.

## What’s next

- **A4 · Choosing a grid recipe** — the four discrete-grid recipes that
  govern the support of $`g`$ and when each one is appropriate. The
  default `paper_realdata` recipe is the one a ten-site analysis like
  this vignette would adopt; the oracle-aware recipes are for simulation
  work where the true $`g`$ is known.
- **A5 · Diagnostics in practice** — the release-quality checklist for
  the
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  output, including the conventional $`\widehat{R} < 1.01`$ and bulk-ESS
  $`> 400`$ thresholds and the remediation patterns for divergences and
  treedepth hits.
- **A7 · Case study — Lee–Sui replication light** — the long-form case
  study that picks up where this vignette leaves off, taking a
  substantive dataset through the full pre-registered analysis and the
  methods-appendix template.
