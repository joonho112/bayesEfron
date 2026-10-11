# Plot sampler diagnostics and grid diagnostics

`plot_diagnostics()` plots the diagnostics returned by
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md):
R-hat and the bulk and tail effective sample sizes (ESS) of the
monitored quantities against their thresholds, and two diagnostics of
the grid. The caption gives the number of divergent transitions, the
number of iterations that reached the maximum tree depth, and the
smallest E-BFMI.

## Usage

``` r
plot_diagnostics(x, what = c("rhat", "ess", "grid"), by_block = TRUE, ...)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- what:

  Character vector with the diagnostics to plot: one or more of
  `"rhat"`, `"ess"` and `"grid"`. By default all three are plotted.

- by_block:

  If `TRUE` (the default), the quantities are grouped into blocks as
  described in Details. If `FALSE`, they are shown in a single row.

- ...:

  These dots are for future extensions and must be empty.

## Value

A `ggplot` object.

## Details

By default the monitored quantities are grouped into blocks: the spline
coefficients `alpha`, the precision `lambda`, the mean, variance and
standard deviation of \\g\\, `log_g` at the grid points, the site-level
quantities (if the fit stores them, see `site_draws` in
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)),
and all other quantities. A block of up to 30 quantities is shown as
points. A larger block is shown as a box plot, with its worst value as a
point. A name next to a block identifies its worst quantity, i.e., the
one with the largest R-hat or the smallest effective sample size.

Dashed lines mark the thresholds used by
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md):
1.01 for R-hat and 400 for the effective sample sizes. A dotted line
marks an R-hat of 1.05, above which the draws should not be used.

Sampler diagnostics are computed from all iterations of the run, even if
the fit was thinned (`thin`) and stores only some of them. The two grid
panels are computed from the stored draws.

Sampler checks that failed or were skipped are listed in the caption.
For a short run (fewer than 400 draws),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
skips the checks of R-hat and of the effective sample sizes. Their
values are still plotted against the thresholds, but a value beyond a
threshold then does not mean that a check failed. Values that are not
available are not plotted, and the caption names the quantities
concerned (a long list is abbreviated). This includes a quantity with an
infinite R-hat: its value is stored as `NA`, and the caption also
reports it as a failed check.

If `what` includes `"grid"`, two more panels show the grid diagnostics.
The first gives, for each end of the grid (its lower and upper limit),
the largest posterior mass that any site has in the outer 2% of the
range, and names that site. A value above the dashed line (0.05) means
that the site has substantial posterior mass at that limit. The second
panel compares the posterior interval of the standard deviation of \\g\\
with its prior interval on the same grid. See the Details of
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
for how to read both panels.

## See also

[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md),
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)

## Examples

``` r
plot_diagnostics(krw_fit)


# R-hat alone, without blocks.
plot_diagnostics(twin_towers_fit, what = "rhat", by_block = FALSE)

```
