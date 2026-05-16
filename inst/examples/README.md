# Cached examples

This directory ships fixtures and scripts used by example blocks in
`man/*.Rd`, vignettes, and the package README.

## `cached_fit_re_smoke.rds`

A small `bef_fit_re` object produced by `bayes_efron_fit()` against a
five-site toy dataset. The fixture is used by the `@examples` blocks on
`bayesEfron-methods`, `diagnose`, and `plot.bef_fit_re` to demonstrate
the S3 surface without requiring `R CMD check` to invoke CmdStan.

**Configuration:**

| Parameter | Value |
|:----------|:------|
| `theta_hat` | `c(-0.21, 0.04, 0.19, 0.38, 0.61)` |
| `sigma` | `c( 0.18, 0.15, 0.22, 0.19, 0.24)` |
| `L` | 51 |
| `M` | 3 |
| `chains` | 1 |
| `iter_warmup` | 150 |
| `iter_sampling` | 100 |
| `seed` | 1234 |
| `grid_method` | `"paper_realdata"` (default) |

**Output dimensions:**

| Field | Size |
|:------|:-----|
| In-memory object size | ~237 KB |
| On-disk `.rds` size | ~140 KB |
| Sites (`K`) | 5 |
| Metadata fields | 13 (per fit-pipeline contract) |

The fixture is intentionally smoke-scale: the small sample size and
short chain produce diagnostics that flag low effective sample size,
which is expected and is the same calibration regime as the package's
release tests. Do not interpret the cached fit as a release-quality
calibration result; it exists to give the S3-method examples something
to operate on.

## `example-bayes-efron-fit.R`

The end-to-end Maya-path example script: starts from group treatment
and control summaries for ten sites, computes effect sizes with
`metafor::escalc()`, converts to a `bef_data` object, fits the model
under a live CmdStan gate, and produces summary / interval /
diagnostic / plot outputs. Sourced by the `bayesEfron-metafor.Rmd`
vignette.

## Regenerating the cached fixture

If `R/`, `inst/stan/`, or grid construction changes in a way that
changes the posterior, regenerate the cached fixture from this
directory with:

```r
# From the package root, with CmdStan installed and configured.
BAYESEFRON_RUN_LIVE=1 Rscript -e '
  pkgload::load_all()
  fit <- bayes_efron_fit(
    theta_hat     = c(-0.21, 0.04, 0.19, 0.38, 0.61),
    sigma         = c( 0.18, 0.15, 0.22, 0.19, 0.24),
    L             = 51L,
    M             = 3L,
    chains        = 1L,
    iter_warmup   = 150L,
    iter_sampling = 100L,
    seed          = 1234L
  )
  saveRDS(fit, "inst/examples/cached_fit_re_smoke.rds")
'
```

The seed is fixed (`1234L`) so the regenerated fixture is
reproducible under the same `R` / CmdStan / package versions.
