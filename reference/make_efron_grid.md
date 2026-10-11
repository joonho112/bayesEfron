# Build the grid and spline basis for the distribution of effects

The model represents the distribution of true effects, \\g\\, by its
probabilities on a grid of equally spaced points. The logarithm of these
probabilities is, up to a normalizing constant, a linear combination of
natural cubic spline functions (see
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)).
`make_efron_grid()` builds the grid and the spline basis from the
estimates.
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
builds its grid in the same way from its own arguments, so you only need
to call this function yourself to inspect a grid before fitting, or to
use one grid in several fits by passing it as `grid`.

## Usage

``` r
make_efron_grid(
  theta_hat,
  sigma = NULL,
  ...,
  L = NULL,
  expansion = 0.5,
  M = 6L,
  grid_method = "pilot",
  theta_true = NULL,
  bound_expansion = NULL,
  support = NULL,
  component = c("point", "normal"),
  width = 2/3,
  knots = "grid",
  margin = 0.1,
  seed = NULL
)

# S3 method for class 'bef_grid'
format(x, digits = 3L, ...)

# S3 method for class 'bef_grid'
print(x, ...)
```

## Arguments

- theta_hat:

  Numeric vector of effect estimates, one per site. Can also be a data
  frame or a `bef_data` object that contains the estimates and their
  standard errors (see
  [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)),
  in which case `sigma` is omitted.

- sigma:

  Numeric vector of standard errors of the estimates. All must be
  positive.

- ...:

  For `make_efron_grid()`, these dots are for future extensions and must
  be empty; the arguments after them must be named. The
  [`print()`](https://rdrr.io/r/base/print.html) method passes them on
  to [`format()`](https://rdrr.io/r/base/format.html).

- L:

  Number of grid points, an integer from 51 to 401. If `NULL` (the
  default), the number is chosen from the data: 101, or more if 101
  points are too few for the standard errors or for the spread of the
  effects. See Details.

- expansion:

  Expansion factor, a number between 0 and 5. With
  `grid_method = "pilot"` or `"paper_realdata"`, the range of
  `theta_hat` is extended on each side by `expansion` times its width.

- M:

  Number of spline functions, an integer from 3 to 10. More spline
  functions allow \\g\\ to have narrower features. In simulated data
  sets of 1,500 sites, the site intervals contained somewhat fewer of
  the true effects with 8 functions than with 6 (see
  [`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)).

- grid_method:

  Method for choosing the limits of the grid: `"pilot"` (the default),
  `"paper_realdata"`, `"paper_simulation"` or `"paper_sensitivity"`. See
  Details. `"pilot"` requires CmdStan; the others do not.

- theta_true:

  Numeric vector of true effects, the same length as `theta_hat`.
  Required by `"paper_simulation"` and `"paper_sensitivity"`, which are
  intended for simulations.

- bound_expansion:

  Expansion factor for `grid_method = "paper_sensitivity"`, a number
  greater than 0 and at most 5: the range of `theta_true` is extended on
  each side by `bound_expansion` times its width. `NULL` (the default)
  is equivalent to 0.5.

- support:

  Numeric vector of length 2 with the lower and upper limits of the
  grid. If `NULL` (the default), the limits are chosen by `grid_method`.

- component:

  Type of component at each grid point: `"point"` (the default) for a
  point mass, or `"normal"` for a normal distribution.

- width:

  Standard deviation of the normal components in units of the grid
  spacing, a number greater than 0 and at most 2. Only used with
  `component = "normal"`.

- knots:

  Method for placing the knots of the spline functions. Currently only
  `"grid"` is available; see Details.

- margin:

  Margin for `grid_method = "pilot"`, a number between 0 and 1. The
  interval that contains the central 99.9% of the pilot estimate of
  \\g\\ is extended on each side by `margin` times its width.

- seed:

  Integer seed for the pilot fit, only used with
  `grid_method = "pilot"`. If `NULL` (the default), a seed is taken from
  the clock.

- x:

  A `bef_grid` object created by `make_efron_grid()`.

- digits:

  Number of significant digits to print.

## Value

An object of class `bef_grid`, which is a list with the following
components.

- `grid`:

  The grid points.

- `B`:

  The spline basis, a matrix with one row per grid point and `M`
  columns.

- `L`, `M`:

  The number of grid points and the number of spline functions.

- `support`:

  The lower and upper limits of the grid.

- `spacing`:

  The distance between neighboring grid points.

- `grid_method`:

  The method that chose the limits, or `"support"` if they were
  specified directly.

- `expansion`:

  The factor by which the range was extended (`NA` if the method does
  not use one).

- `component`:

  `"point"` or `"normal"`.

- `width`:

  The value of `width` (`NA` for point masses).

- `tau_c`:

  The standard deviation of a component (0 for point masses).

- `knots`:

  A list with the method and the positions of the interior and boundary
  knots.

- `resolution`:

  A list that compares the grid spacing with the data, as described in
  Details. It contains the two ratios, `sigma_min_over_spacing` and
  `scale_over_spacing`; `L_needed`, the smallest number of points for
  which both are large enough; `adequate`, which is `TRUE` if they are
  large enough for this grid; `L_chosen`, which is `TRUE` if the number
  of points was chosen from the data; `share_outside_estimates`, the
  proportion of grid points outside the range of the estimates; and
  `width_over_sigma`, the width of the grid divided by the median
  standard error.

- `attribution`:

  A list with the formula of the method (`formula`) and its source
  (`source`).

- `pilot`:

  Only for `grid_method = "pilot"`: a list that records what the method
  did. `wide` is the initial range and `core` is the interval that
  contained the central 99.9% of the pilot estimate of \\g\\. `narrowed`
  is `TRUE` if the range was narrowed, `share` is the proportion of the
  initial range that remains, and `protected` contains the sites whose
  estimates caused a limit to be moved outward. `margin` and `seed` are
  the values that were used. `chains`, `iter_warmup` and `iter_sampling`
  give the length of the pilot fit (zero if none was run), `diagnostics`
  contains its sampler diagnostics, and `seconds` is its run time.

- `data`:

  A list with the estimates and standard errors from which the grid was
  built, used by
  [`plot.bef_grid()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid.md).
  It is not included in the grid stored in a fit.

## Details

The grid is determined by four choices.

- The range, which is either chosen by the method named in `grid_method`
  or specified directly as `support`.

- The number of grid points, `L`.

- What \\g\\ places at each grid point, `component`. With `"point"` (the
  default, and the model of Lee and Sui, 2025) it is a point mass. With
  `"normal"` it is a normal distribution whose standard deviation is
  `width` times the grid spacing, as in the penalized Gaussian mixture
  of Sui and Tipton (2026).

- The flexibility of \\g\\, which depends on the number of spline
  functions, `M`, and on the placement of their knots, `knots`.

There are four methods for choosing the range. Three of them come from
the code that accompanies Lee and Sui (2025), and `"pilot"` narrows the
range of the first of them.

|  |  |  |
|----|----|----|
| `grid_method` | Limits of the grid | Needs `theta_true` |
| `"paper_realdata"` | The range of `theta_hat`, extended on each side by `expansion` times its width. | No |
| `"pilot"` | The range of `"paper_realdata"`, narrowed to the region where a short preliminary fit on that range (the pilot fit) places \\g\\. | No |
| `"paper_simulation"` | The range of `theta_true`, extended by 0.5 (in the units of the effects) on each side. | Yes |
| `"paper_sensitivity"` | The range of `theta_true`, extended on each side by `bound_expansion` times its width. | Yes |

Use `"pilot"` (the default) or `"paper_realdata"` for real data;
`"paper_realdata"` does not require a model fit. The other two methods
need the true effects, so they can only be used in simulations. They are
provided so that the paper's simulations can be reproduced.

With the default `expansion`, the range of `"paper_realdata"` is twice
as wide as the range of the estimates, and the estimates are themselves
more spread out than the true effects. With many sites, \\g\\ then
occupies only a small part of the grid, most of the spline functions
describe the empty part, and sampling is inefficient. The `"pilot"`
method avoids this in two steps.

1.  A pilot fit is run on the range of `"paper_realdata"`, with 4 chains
    of 500 warmup and 500 sampling iterations, 101 grid points and 6
    spline functions. The pilot fit always uses point masses, even when
    the final grid has normal components.

2.  The method finds the grid points that enclose the central 99.9% of
    the pilot estimate of \\g\\ and extends this interval on each side
    by `margin` times its width, and by at least one grid spacing of the
    pilot fit. If an estimate lies more than 4 of its standard errors
    outside a limit, the limit is moved out to that estimate plus or
    minus 2 standard errors.

The limits never go beyond the range of `"paper_realdata"`, so the
method either narrows that range or leaves it unchanged. It leaves the
range unchanged for
[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
and
[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md)
and narrows it for
[twin_towers](https://joonho112.github.io/bayesEfron/reference/twin_towers.md);
see
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md).
If the range of `"paper_realdata"` is narrower than the median standard
error, no pilot fit is run. The method requires CmdStan. Specify `seed`
to make the result reproducible. The `pilot` component of the result
records what the method did.

If `support = c(lower, upper)` is specified, the grid runs from `lower`
to `upper` and no method is applied, so `grid_method` must not be
specified. With point masses, \\g\\ has no probability outside the
support, so the support is the set of values that an effect can take in
the model. With normal components, the support is the range of the
component centers, and \\g\\ extends about 2 component standard
deviations further on each side. If the effects must stay within given
bounds, use point masses together with `reading = "step"` in the
summaries, because with the default `reading` a summary can extend half
a grid spacing beyond a limit (see
[`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md)).

With `knots = "grid"`, the spline basis is
`splines::ns(grid, df = M, intercept = FALSE)`. The interior knots are
at quantiles of the grid points, which makes them equally spaced over
the grid, and the boundary knots are at its limits. An error is thrown
if the basis together with a constant column does not have full column
rank; a smaller `M` avoids this.

The grid must also be fine enough for the data. With point masses, the
likelihood of a site is evaluated only at the grid points, so a site
whose standard error is much smaller than the grid spacing is
represented poorly. In addition, a grid with only a few points across
the bulk of the effects cannot show the shape of \\g\\. The function
computes two ratios to check this.

- The smallest standard error divided by the grid spacing, which should
  be at least 0.75.

- A rough estimate of the standard deviation of the effects divided by
  the grid spacing, which should be at least 4. With \\m\\ the median
  absolute deviation of the estimates (scaled as in
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html), so that it
  estimates a standard deviation) and \\v\\ the median sampling
  variance, this estimate is \\\sqrt{m^2 - v}\\, or \\\sqrt{v}\\ if
  \\m^2 - v\\ is not positive.

With `L = NULL`, the grid has 101 points if both ratios are large
enough, and otherwise the smallest number of points, up to 401, for
which they are. With normal components only the second ratio is used, in
an equivalent form: the standard deviation of a component should be at
most one sixth of the estimated standard deviation of the effects. A
standard error smaller than the grid spacing is not a problem in that
case. A warning of class `bef_warn_grid_resolution` is issued if a ratio
is too small for the returned grid, or if the whole grid is narrower
than the median standard error. The warning message suggests what to
change. The ratios are stored in the `resolution` component of the
result.

## Functions

- `format(bef_grid)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows as a character
  vector.

- `print(bef_grid)`: Prints the limits and the spacing of the grid, the
  method that chose the limits, the components, the spline functions
  with their knots, and how the spacing compares with the data.

## References

Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
2639. [doi:10.3390/math13162639](https://doi.org/10.3390/math13162639)

Sui, D. and Tipton, E. (2026). Meta-analysis without normality:
estimating the true effect distribution with penalized Gaussian
mixtures. arXiv:2604.27887.

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md),
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)

## Examples

``` r
# The default method runs the pilot fit, which needs CmdStan.
if (FALSE) { # \dontrun{
make_efron_grid(twin_towers$yi, twin_towers$sei, seed = 20250627)
} # }

# The method "paper_realdata" needs no fit.
grid <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata")
grid
#> <bef_grid> 101 points from -0.0829 to 0.159, spacing 0.00241
#>   Rule:       "paper_realdata", expansion 0.5
#>   Component:  point masses
#>   Basis:      6 natural cubic spline functions, knots by rule "grid"
#>   Interior knots: -0.0427, -0.00242, 0.0378, 0.0781, 0.118
#>   Resolution: smallest standard error 2.21 spacings, scale of the effects 7.89 spacings
#>   50.5% of the grid points lie outside the range of the estimates
range(krw_discrimination$yi)
#> [1] -0.02254064  0.09817050
grid$support
#> [1] -0.08289621  0.15852607
dim(grid$B)
#> [1] 101   6

# The six spline functions and the knots over the estimates, in
# percentage points.
plot(grid, scale = 100)


# A wider grid with more spline functions.
wide <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata",
                        expansion = 1, M = 8)
wide$support
#> [1] -0.1432518  0.2188816

# The limits specified directly, and a normal distribution at each grid
# point.
make_efron_grid(krw_discrimination, support = c(-0.05, 0.15),
                component = "normal")
#> <bef_grid> 101 points from -0.05 to 0.15, spacing 0.002
#>   Rule:       ends given in `support`
#>   Component:  normal, standard deviation 0.00133 (0.667 of the spacing)
#>   Basis:      6 natural cubic spline functions, knots by rule "grid"
#>   Interior knots: -0.0167, 0.0167, 0.05, 0.0833, 0.117
#>   Resolution: smallest standard error 2.67 spacings, scale of the effects 9.52 spacings
#>   39.6% of the grid points lie outside the range of the estimates

# In a simulation the true effects can set the range.
sim <- make_efron_grid(twin_towers, grid_method = "paper_simulation",
                       theta_true = twin_towers$theta_true)
sim$support
#> [1] -2.198691  3.197622
```
