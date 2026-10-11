# Plots of a fitted model

This vignette shows what each plot function in bayesEfron is for and how
to read the result. A fit contains a posterior distribution for \\g\\,
the distribution of effects, and for the effect of every site, and most
questions about it are easier to answer from a plot than from a table.

``` r

library(bayesEfron)
```

Each plot function returns a `ggplot` object that you can change with
ggplot2 functions. The one exception is `plot_mcmc(type = "pairs")`,
which returns the grid of plots that bayesplot makes.

For most examples we use two fits stored in the package: `krw_fit`, the
callback gaps of 97 firms from
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md),
and `twin_towers_fit`, 1,500 simulated sites with known true effects
from
[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md).

``` r

firms <- krw_fit
sims <- twin_towers_fit
truth <- twin_towers$theta_true
```

## The distribution of effects

[`plot_g()`](https://joonho112.github.io/bayesEfron/reference/plot_g.md)
shows what the distribution of true effects looks like and how well the
data determine it. The curve is the posterior mean of the density of
\\g\\, and the pointwise band around it contains the density with the
stated probability at each value. Below, `observed` puts a histogram of
the estimates under the curve, `normal` puts a normal distribution
fitted to the same estimates over it, and `scale = 100` gives percentage
points.

``` r

plot_g(firms, observed = "histogram", normal = TRUE, scale = 100)
```

![The estimated density of the callback gaps with its band, a histogram
of the estimates and a dashed normal
curve.](plots_files/figure-html/g-1.png)

The histogram is wider than the curve because the estimates contain
sampling error and the effects do not. Where the band is wide, the data
do not determine the height of the density well.

Proportions of sites are easier to read from the distribution function,
so next we set `type = "cdf"`. Setting `inner_level` adds a second,
narrower band.

``` r

plot_g(firms, type = "cdf", inner_level = 0.5, scale = 100)
```

![The estimated distribution function of the callback gaps with a 50
percent band inside a 90 percent
band.](plots_files/figure-html/g-cdf-1.png)

The height of the curve at a value is the proportion of firms with a gap
below it. Three more arguments add to the plot:

- `draws = 40` adds 40 posterior draws of \\g\\ to its posterior mean
  and band.
- `truth` adds the true effects of simulated data.
- `compare` puts other estimates of the distribution over this one, from
  other fits or from other software.

## The sites

### Caterpillar plot

A caterpillar plot from
[`plot_sites()`](https://joonho112.github.io/bayesEfron/reference/plot_sites.md)
shows the posterior mean of each site with a 50% and a 90% interval, and
the estimate as an open circle. Use it to see which sites have the
largest and the smallest effects, and how certain that ordering is.

``` r

plot_sites(firms, n_extreme = 10, scale = 100)
```

![Two panels with the ten firms with the highest and the ten with the
lowest posterior mean gaps, by
name.](plots_files/figure-html/sites-1.png)

`n_extreme` selects the sites at the two ends, and `sites` selects sites
by name or position. Without either, every site is plotted: labeled when
there are at most 100, and as thin lines when there are more.

``` r

plot_sites(sims)
```

![The 90 percent intervals of the 1,500 simulated sites as thin lines in
the order of the posterior
means.](plots_files/figure-html/sites-many-1.png)

### Shrinkage plot

[`plot_shrinkage()`](https://joonho112.github.io/bayesEfron/reference/plot_shrinkage.md)
shows how far the model moves the estimates by plotting the posterior
mean of every site against its estimate. Points on the diagonal have not
moved.

``` r

plot_shrinkage(firms, scale = 100)
```

![The posterior mean of each firm against its estimate, colored by
standard error.](plots_files/figure-html/shrinkage-1.png)

With `style = "arrows"`, each site is an arrow from the estimate to the
posterior mean, at the height of its standard error. This suits a small
number of sites, such as the 19 teacher expectancy studies in the stored
fit `raudenbush_fit`.

``` r

plot_shrinkage(raudenbush_fit, style = "arrows")
```

![For each of the 19 teacher expectancy studies an arrow from the
estimate to the posterior mean, at the height of its standard
error.](plots_files/figure-html/shrinkage-arrows-1.png)

### Exceedance probabilities

When the question is which sites have an effect above a threshold,
[`plot_exceedance()`](https://joonho112.github.io/bayesEfron/reference/plot_exceedance.md)
plots, for each site, the posterior probability that its effect is above
the threshold. `cutoff` marks a probability and counts the sites that
reach it. For many sites, the probabilities are plotted as a curve.

``` r

plot_exceedance(sims, cutoff = 0.9)
```

![The probabilities that the effects of the 1,500 simulated sites are
positive, in increasing order, with a dashed line at
0.9.](plots_files/figure-html/exceedance-1.png)

The caption notes that the count controls no error rate over the sites,
which is easy to forget when a list of sites is made from it.

### The posterior of one site

[`plot_site_posterior()`](https://joonho112.github.io/bayesEfron/reference/plot_site_posterior.md)
plots the whole posterior distribution of a site, with the two things
that it is made from: the likelihood of the estimate and the estimate of
\\g\\.

``` r

plot_site_posterior(sims, sites = c(268, 304))
```

![Two panels: an imprecise site between the two groups, whose posterior
distribution has two modes, and the most precise site, whose posterior
distribution is narrow.](plots_files/figure-html/site-posterior-1.png)

The first site has an imprecise estimate between the two groups of
simulated effects, and its posterior distribution has two modes. A mean
and an interval would hide that.

## Checking a fit

### Sampler and grid diagnostics

[`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md)
plots what
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
reports, so you can check two things at once: whether the sampler
converged, and whether the grid suits the fit. R-hat and the effective
sample sizes of the quantities that the sampler records are shown in
blocks against their thresholds, with two more panels about the grid.

``` r

plot_diagnostics(firms)
```

![Five panels: R-hat, bulk and tail effective sample sizes by block, the
largest mass of a site at each end of the grid, and the standard
deviation of the distribution of effects under the prior and the
posterior.](plots_files/figure-html/diagnostics-1.png)

Points to the right of the dashed line in the first panel (R-hat), or to
the left of it in the next two (effective sample sizes), would be
quantities that miss their thresholds. There are none here. The caption
gives the numbers of divergent transitions and of iterations at the
maximum tree depth, and
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
explains each check.

### Posterior predictive check

A posterior predictive check with
[`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md)
shows whether the model reproduces the estimates. The check compares
each estimate with the distribution that the fitted model gives it,
through the probability of a smaller estimate at that site (the
probability integral transform, PIT). If \\g\\ were known, these values
would be uniformly distributed between 0 and 1, and the line shows how
far the values of the fitted model are from that.

``` r

plot_ppc(firms)
```

![The empirical distribution function of the probability integral
transforms of the 97 estimates less the uniform one, inside a gray
band.](plots_files/figure-html/ppc-1.png)

The band is the reference for independent uniform values, but it does
not allow for the fitting of \\g\\. See
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
for how to read it.

### Pareto k values

Leave-one-out cross-validation approximates the fit that would result
without each site, and the Pareto k of a site measures how good that
approximation is. To check that the approximation is reliable for every
site, plot the values with
[`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md).
A line marks the value above which the result for a site should not be
trusted, and the sites above it are labeled. See
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
for more.

``` r

plot_loo(firms)
```

![The Pareto k value of each of the 97 firms, all below the dashed
line.](plots_files/figure-html/loo-1.png)

## Grids and grid comparisons

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
a grid, or `plot(fit, type = "grid")` for the grid of a fit, shows where
the grid lies in relation to the data: the spline functions, the knots
and the grid points over the estimates.

``` r

plot(sims, type = "grid")
```

![The eight spline basis functions and the knots of the grid of the
simulated-site fit, over a histogram of the
estimates.](plots_files/figure-html/grid-1.png)

The limits of this grid are well inside the range of the estimates
because the default `grid_method` narrowed it.

To see what changes when the grid changes, plot a comparison from
[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md):
the specifications appear one under the other. Here we load a stored
comparison of three ranges for the 19 teacher expectancy studies.

``` r

comparison <- readRDS(system.file("examples", "raudenbush_ranges.rds",
                                  package = "bayesEfron"))
plot(comparison)
```

![Three rows, one for each range of the grid for the 19 studies, each
with the estimated density, the mean and standard deviation of the
distribution, and the posterior means and intervals of two
studies.](plots_files/figure-html/comparison-1.png)

[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)
discusses the first of these plots and
[`vignette("teacher-expectancy")`](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md)
the second.

## Simulated data

With simulated data,
[`plot_recovery()`](https://joonho112.github.io/bayesEfron/reference/plot_recovery.md)
shows how close a fit comes to the known effects. It plots the intervals
of the sites with the true effects and marks the sites where the true
effect is outside. Only 150 sites are shown, and `seed` fixes the random
choice.

``` r

plot_recovery(sims, truth, seed = 1)
```

![The 90 percent posterior intervals of 150 sites with their true
effects as crosses, and as orange triangles where they lie outside the
interval.](plots_files/figure-html/recovery-1.png)

Setting `type = "calibration"` shows the proportion of true effects
inside the intervals at every level from 5% to 95%, and
`type = "scatter"` plots the posterior means against the true effects.

## Common arguments

Every plot function takes the fit first, and a few arguments work the
same way wherever they apply:

- `scale` multiplies the effect values. For example, `scale = 100` turns
  differences in probabilities into percentage points.
- `level` sets the probability of posterior intervals and bands. In
  [`plot_ppc()`](https://joonho112.github.io/bayesEfron/reference/plot_ppc.md)
  it sets the level of the independent-uniform reference intervals.
- `sites` and `highlight` take positions or labels.

Calling [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
fit with a `type` runs the same functions by name, so
`plot(firms, type = "g", scale = 100)` is the same as
`plot_g(firms, scale = 100)`.

## Customizing a plot

A plot has no title, but its caption says what the lines and bands are.
You can set both, and add layers, as with any other `ggplot` object:

``` r

plot_g(firms, scale = 100) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
  ggplot2::labs(title = "Callback gaps of 97 firms",
                x = "White minus Black callback rate (percentage points)")
```

![The density of the callback gaps with a title, an axis label and a
vertical line at zero added.](plots_files/figure-html/layers-1.png)

## Making your own plots

The numbers in the plots of \\g\\ and of the sites are also available
from three functions. The densities, distribution functions and
probabilities of \\g\\ come from
[`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md),
the means and intervals of the sites from
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
and the probabilities of exceeding a threshold from
[`exceedance_prob()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md).
To make a plot that the package does not have, start from these
functions:

``` r

head(g_summary(firms, type = "cdf", scale = 100), 3)
#>       value         mean        lower        upper
#> 1 -8.289621 0.0002198804 2.447210e-07 0.0008876288
#> 2 -8.048199 0.0005935610 7.766047e-07 0.0024197491
#> 3 -7.806777 0.0008606321 1.298678e-06 0.0034986411
site_summary(firms, sites = c("KFC", "Walmart"), scale = 100)[, 1:8]
#>   site   label estimate       se     mean       sd      lower    upper
#> 1    2     KFC      2.2 1.530472 1.786494 1.052594 0.12271869 3.581179
#> 2   65 Walmart      7.0 3.919826 2.386823 1.744476 0.04493205 5.523267
```

## MCMC plots with bayesplot

Trace plots, rank plots and the other plots of the Markov chains come
from bayesplot, and
[`plot_mcmc()`](https://joonho112.github.io/bayesEfron/reference/plot_mcmc.md)
passes the draws of a fit to it. We show a rank plot, in which the draws
of all the chains are ranked together and each line counts the ranks of
one chain in 20 bins. Chains that sample the same distribution give
lines that overlap around a common level, as they do here.

``` r

plot_mcmc(firms, type = "rank", pars = c("lambda", "sd_g"))
```

![Rank plots of lambda and of the standard deviation of the distribution
of effects: for each, four lines, one for each chain, that overlap
around a common level.](plots_files/figure-html/mcmc-1.png)

The `type` argument can be `"trace"`, `"rank"`, `"dens"`, `"pairs"` or
`"energy"`, and
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
shows a trace plot. For anything else, pass the draws and the sampler’s
quantities to bayesplot yourself:

``` r

bayesplot::mcmc_areas(firms$draws, pars = c("mean_g", "sd_g"))
bayesplot::mcmc_nuts_divergence(bayesplot::nuts_params(firms),
                                bayesplot::log_posterior(firms))
```

Note that the stored fits are thinned to 1 draw in 15, so their trace
plots show 200 draws for each chain, whereas the sampler diagnostics in
[`plot_diagnostics()`](https://joonho112.github.io/bayesEfron/reference/plot_diagnostics.md)
are computed from the whole run. The energy plot compares each iteration
with the one before it and is available only for a fit that stores every
iteration. See
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
for the call.
