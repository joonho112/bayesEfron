# Convert input data to a bayesEfron data object

Adapt the user's effect-size data into the standalone `bef_data` class
that
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
consumes. Three input shapes are supported, all returning a validated
`bef_data` object that carries `theta_hat`, `sigma`, optional site
labels, and a `source` attribute recording the input shape that produced
it.

## Usage

``` r
as_bef_data(x, ...)

# Default S3 method
as_bef_data(x, ...)

# S3 method for class 'list'
as_bef_data(x, ...)

# S3 method for class 'escalc'
as_bef_data(x, ...)
```

## Arguments

- x:

  Object to convert. v0.1 supports `bef_data`, `list` with `theta_hat`
  and `sigma`, and
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  objects.

- ...:

  Reserved for future expansion; must be empty in v0.1.

## Value

A validated `bef_data` object with a `source` attribute set to one of
`"list"`, `"metafor::escalc"`, or unchanged for already converted
inputs.

## Details

Supported input shapes:

|  |  |  |
|----|----|----|
| Input class | Required fields | Source label |
| `bef_data` | (already converted; revalidated and returned) | unchanged |
| named `list` | `theta_hat`, `sigma`; optional `names` | `"list"` |
| `escalc` from [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html) | `yi`, `vi` (variance); optional row labels | `"metafor::escalc"` |

For the `escalc` path, `theta_hat` is taken from `yi` and `sigma` from
`sqrt(vi)` after a strict-positivity check on `vi`. Row names on the
`escalc` object are propagated as site labels when present and
non-default.

For the `list` path, an explicit `names` element overrides any names
attribute on `theta_hat`. Site labels are dropped when they are missing,
empty, or contain `NA`.

Unsupported input classes raise a typed `bef_error` via
`as_bef_data.default()` with the offending class name in the message.

## See also

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
  which calls `as_bef_data()` internally on its `theta_hat` / `sigma`
  arguments.

- [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  for computing effect sizes and sampling variances from study-level
  summary data.

## Examples

``` r
# Plain list input.
dat_list <- as_bef_data(list(
  theta_hat = c(-0.21, 0.04, 0.19, 0.38, 0.61),
  sigma     = c( 0.18, 0.15, 0.22, 0.19, 0.24)
))
dat_list
#> <bef_data>
#> Sites: 5
#> Source: list
#> theta_hat: min -0.21; median 0.19; max 0.61
#> sigma: min 0.15; median 0.19; max 0.24

# Optional site labels (must satisfy the minimum-length-5 constraint).
dat_named <- as_bef_data(list(
  theta_hat = c(
    site_1 = -0.21, site_2 = 0.04, site_3 = 0.19,
    site_4 =  0.38, site_5 = 0.61
  ),
  sigma = c(0.18, 0.15, 0.22, 0.19, 0.24)
))

# metafor::escalc() bridge.
if (requireNamespace("metafor", quietly = TRUE)) {
  esc <- metafor::escalc(
    measure = "MD",
    m1i  = c(0.10, 0.40, 0.55, 0.30, 0.20),
    sd1i = c(0.30, 0.30, 0.35, 0.28, 0.32),
    n1i  = c( 60,   55,   62,   58,   65),
    m2i  = c(0.05, 0.10, 0.15, 0.08, 0.12),
    sd2i = c(0.32, 0.34, 0.36, 0.30, 0.33),
    n2i  = c( 60,   55,   62,   58,   65)
  )
  dat_esc <- as_bef_data(esc)
  dat_esc
}
#> <bef_data>
#> Sites: 5
#> Source: metafor::escalc
#> theta_hat: min 0.05; median 0.22; max 0.4
#> sigma: min 0.05388; median 0.05702; max 0.06377
```
