# Plot the sensitivity of the results to the grid

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
a comparison created by
[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
shows one row per specification, so that you can see how the results
change with the grid. Each row shows the estimated density of \\g\\ with
its pointwise posterior band, the posterior means of the mean and of the
standard deviation of \\g\\, and the posterior means and intervals of
the effects of the selected sites. Read down a column to see what
changes when the range, the number of grid points or the number of
spline functions is varied.

## Usage

``` r
# S3 method for class 'bef_grid_comparison'
plot(x, scale = 1, ...)
```

## Arguments

- x:

  A `bef_grid_comparison` object created by
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md).

- scale:

  Multiplier for the values on the effect axis, a positive number. For
  example, `scale = 100` shows differences in proportions as percentage
  points.

- ...:

  These dots are for future extensions and must be empty.

## Value

A `ggplot` object.

## Details

The specifications are shown in the order in which they were fitted,
with the first at the top. They are labeled by the settings that differ
between them; in these labels, "rule" means the grid method. A
specification whose fit failed still has a row, which is labeled as
failed and is otherwise empty.

All densities are plotted on the same vertical scale, so their heights
can be compared between rows, but the scale itself is not shown. To
avoid exaggerating small differences between specifications, the axis
for the mean of \\g\\ is at least as wide as the standard deviation of
\\g\\, and the axis for the standard deviation starts at zero.

The sites are those selected with the `sites` argument of
[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md).
By default these are the sites with the smallest and the largest
estimates, whose effects depend most on the limits of the grid.

## See also

[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md),
[`plot.bef_grid()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_grid.md),
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)

## Examples

``` r
# Three ranges for the 19 teacher expectancy studies, stored with the
# package. The examples of compare_efron_grids() show the call that made
# the comparison.
comparison <- readRDS(system.file("examples", "raudenbush_ranges.rds",
                                  package = "bayesEfron"))
comparison[, c("expansion", "lower", "upper", "sd_g", "site_4")]
#> <bef_grid_comparison> 3 specifications, best first
#> 
#>  expansion  lower upper  sd_g site_4
#>       0.50 -1.070  1.93 0.314  0.663
#>       0.25 -0.695  1.56 0.317  0.782
#>       1.00 -1.820  2.68 0.297  0.442

# As the range widens, the mean of study 4 moves toward the others.
plot(comparison)

```
