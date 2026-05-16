# Construct an Efron log-spline grid

Build the discrete support and natural-cubic-spline basis that the v0.1
random-effects Stan model uses to represent the mixing distribution
\\g\\. The function exposes the four grid recipes committed by the
package blueprint and the associated MDPI *Mathematics* paper, together
with one experimental KL-target recipe.

Most users do not need to call `make_efron_grid()` directly:
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
constructs a grid internally with the same arguments. Calling this
function on its own is useful when you want to inspect or visualize a
candidate grid before committing to a fit, or to share the same grid
across several fits for comparison.

## Usage

``` r
make_efron_grid(
  theta_hat,
  sigma,
  L = NULL,
  expansion = 0.5,
  kappa = NULL,
  M = 6L,
  grid_method = "paper_realdata",
  theta_true = NULL,
  bound_expansion = NULL
)
```

## Arguments

- theta_hat:

  Numeric vector of observed effect estimates, length \\\ge 2\\.

- sigma:

  Numeric vector of strictly positive within-study standard errors. Same
  length as `theta_hat`.

- L:

  Integer grid length, or `NULL`. `NULL` uses `101L` for the three paper
  recipes; `"kl_target_experimental"` derives `L` from `kappa`. Bounded
  to \\\[51, 300\]\\ when supplied.

- expansion:

  Numeric, range-relative expansion factor used by `"paper_realdata"`
  and `"kl_target_experimental"`. Defaults to `0.5` (50 percent
  expansion). Range \\\[0, 5\]\\.

- kappa:

  Numeric, KL target used by `"kl_target_experimental"`. `NULL` defaults
  to `1 / length(theta_hat)`. Range \\(0, 1)\\.

- M:

  Integer natural-cubic-spline degrees of freedom. Defaults to `6L`.
  Range \\\[3, 10\]\\.

- grid_method:

  Character grid recipe. One of `"paper_realdata"` (default),
  `"paper_simulation"`, `"paper_sensitivity"`, or
  `"kl_target_experimental"`.

- theta_true:

  Numeric oracle vector required by `"paper_simulation"` and
  `"paper_sensitivity"`; ignored by the other recipes. Same length as
  `theta_hat` when supplied.

- bound_expansion:

  Numeric, oracle-bound expansion factor used only by
  `"paper_sensitivity"`. `NULL` falls back to the recipe default of
  `0.5`. Range \\(0, 5\]\\.

## Value

A named list with the following fields:

- `grid` — numeric vector of length `L`, the discrete support of the
  mixing distribution.

- `B` — natural-cubic-spline basis matrix from
  `splines::ns(grid, df = M, intercept = FALSE)`.

- `M` — integer, the spline degrees of freedom actually used.

- `L` — integer, the grid length actually used.

- `expansion` — the range-relative expansion factor applied
  (recipe-dependent).

- `kappa` — the KL target (only meaningful for
  `"kl_target_experimental"`; `NULL` otherwise).

- `grid_method` — the recipe name.

- `attribution` — a list recording the formula and source lineage of the
  recipe.

## Grid recipes

|  |  |  |
|----|----|----|
| Recipe | Needs `theta_true`? | Use when |
| `"paper_realdata"` (default) | No | Real-data analysis with no oracle. Endpoints come from the observed range of `theta_hat`, expanded by `expansion`. |
| `"paper_simulation"` | Yes | Simulation with known truth, matched to the paper's grid rule. Endpoints come from the oracle range of `theta_true`, padded by an absolute 0.5 on each side. |
| `"paper_sensitivity"` | Yes | Sensitivity sweep that widens the oracle bounds by `bound_expansion`. |
| `"kl_target_experimental"` | No | KL-target tuning. Computes `L` from `kappa`; emits a once-per-session disclaimer because the calibration is experimental for heteroscedastic inputs. |

Across recipes the returned grid is always strictly increasing, covers
the observed range, and supports a length-`M` natural-cubic-spline basis
(`splines::ns(grid, df = M, intercept = FALSE)`).

## Defaults and bounds

`L` defaults to `101L` for the three paper-faithful recipes.
`"kl_target_experimental"` always derives the effective `L` from `kappa`
after validating user inputs; a supplied `L` is range-checked but does
not participate in the recipe's effective grid length. Both `L` (51–300)
and `M` (3–10) are bounded; out-of-range values are rejected at the
input boundary. `expansion` and `bound_expansion` accept values in
\\\[0, 5\]\\; `kappa` accepts values in \\(0, 1)\\.

## Attribution

Every returned grid carries an `attribution` slot recording the formula
and source lineage of the recipe so downstream consumers can audit which
paper rule produced the grid.

## See also

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
  which constructs a grid internally and is the usual entry point.

- The methodological vignette M4 (grid construction and spline basis)
  for the mathematical specification of each recipe.

## Examples

``` r
theta_hat <- c(-0.45, -0.10, 0.20, 0.55, 0.90)
sigma     <- c( 0.20,  0.18, 0.22, 0.16, 0.24)

# Default real-data grid: endpoints from observed range + 50 percent expansion.
g_real <- make_efron_grid(theta_hat, sigma)
length(g_real$grid)            # 101 by default
#> [1] 101
dim(g_real$B)                  # L x M
#> [1] 101   6
g_real$grid_method
#> [1] "paper_realdata"

# Simulation grid (oracle required): pretend the truth is known.
theta_true <- c(-0.40, 0.00, 0.10, 0.70, 1.00)
g_sim <- make_efron_grid(
  theta_hat   = theta_hat,
  sigma       = sigma,
  theta_true  = theta_true,
  grid_method = "paper_simulation"
)

# Sensitivity grid: widen oracle bounds by `bound_expansion`.
g_sens <- make_efron_grid(
  theta_hat       = theta_hat,
  sigma           = sigma,
  theta_true      = theta_true,
  grid_method     = "paper_sensitivity",
  bound_expansion = 0.5
)

# Experimental KL-target recipe: emits a once-per-session disclaimer.
g_kl <- make_efron_grid(
  theta_hat   = theta_hat,
  sigma       = sigma,
  grid_method = "kl_target_experimental"
)
#> kl_target_experimental: ebnm KL bound assumes homoscedastic observations (`grid_selection.R:60-63`); calibration on heteroscedastic input is not validated against any published number.
#> kl_target_experimental grid length set to 51 to enforce [51, 300] bounds.
```
