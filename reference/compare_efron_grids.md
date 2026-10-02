# Compare grids and numbers of spline functions by cross-validation

The number of grid points, the width of the grid and the number of
spline functions are choices of the analyst. `compare_efron_grids()`
fits the model once for every combination of the values given in `L`,
`M` and `expansion` and ranks the fits by `elpd_loo`, the expected log
predictive density for a new site as estimated by leave-one-out
cross-validation
([`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md)).

## Usage

``` r
compare_efron_grids(
  theta_hat,
  sigma,
  L = c(51L, 101L, 201L),
  M = c(4L, 6L, 8L),
  expansion = 0.5,
  grid_method = "paper_realdata",
  ...
)

# S3 method for class 'bef_grid_comparison'
print(x, ...)
```

## Arguments

- theta_hat:

  Numeric vector of effect estimates, one for each site, all on the same
  scale (standardized mean differences or log odds ratios, for example).
  At least five are needed.

- sigma:

  Numeric vector of the standard errors of the estimates, positive and
  of the same length as `theta_hat`.

- L:

  Vector of whole numbers between 51 and 300, the numbers of grid points
  to try. Defaults to `c(51, 101, 201)`.

- M:

  Vector of whole numbers between 3 and 10, the numbers of spline
  functions to try. Defaults to `c(4, 6, 8)`.

- expansion:

  Numeric vector of values between 0 and 5, the expansion factors to
  try. Defaults to `0.5`.

- grid_method:

  Character string, the rule for the end points of the grid:
  `"paper_realdata"` (the default), `"paper_simulation"`,
  `"paper_sensitivity"` or `"kl_target_experimental"`. The Details of
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  describe them.

- ...:

  For `compare_efron_grids()`, other arguments of
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
  such as `chains`, `iter_sampling` and `seed`. Not used by
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  A `bef_grid_comparison` object.

## Value

A data frame of class `bef_grid_comparison` with one row for each
specification, the best first. Its columns are `L`, `M` and `expansion`;
`elpd_loo`, its standard error `se_elpd_loo` and `p_loo`, the effective
number of parameters as the loo package defines it, which is not the
`effective_params` of a fit; `elpd_diff` and `se_diff`; `max_rhat`,
`min_ess_bulk`, `divergences` and `max_pareto_k`; and `sampler_seconds`,
the time that sampling took.

## Details

The table is ordered from the largest `elpd_loo` to the smallest.
`elpd_diff` is the difference from the first row and `se_diff` is its
standard error. A difference of less than about two standard errors does
not separate two specifications. With few sites the standard error is
itself uncertain, and the rule is a rough guide.

A fit that has not converged should not be chosen for its `elpd_loo`.
Each row therefore carries the largest R-hat, the smallest bulk
effective sample size, the number of divergent transitions and the
largest Pareto k of its fit, and
[`print()`](https://rdrr.io/r/base/print.html) names the rows with an
R-hat above 1.05, a divergent transition or a Pareto k above 0.7.

With the default arguments nine models are fitted, each with the sampler
settings of
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
unless they are changed through `...`. A fit that fails gives a warning
and a row of missing values. The function needs the loo package.

## See also

[`loo.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/loo.bef_fit.md),
[`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md),
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)

## Examples

``` r
# Nine fits of the teacher expectancy data. It needs CmdStan.
if (FALSE) { # \dontrun{
compare_efron_grids(raudenbush1985$yi, raudenbush1985$sei, seed = 1985)
} # }
```
