# Reproducibility and interval coverage

If you report results from bayesEfron, you may be asked whether the fit
can be reproduced and what evidence there is that its intervals have
their nominal coverage. Here we collect what you need to answer: how to
reproduce a fit, what a fitted object stores, the two Stan programs, the
checks that we ran on the computations, interval coverage in simulated
data, and how the package relates to the paper ([Lee and Sui
2025](#ref-lee2025)).

``` r

library(bayesEfron)
```

Our example is `krw_fit`, a fit that comes with the package. It is the
fit to the callback gaps of 97 firms that
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)
analyzes.

``` r

fit <- krw_fit
```

## Reproducing a fit

A fitted object records what produced its draws. We first print the
seed, the CmdStan version, the name of the Stan program and the first
characters of a hash of the program’s file.

``` r

fit$metadata$seed
#> [1] 20250807
fit$metadata$cmdstan_version
#> [1] "2.38.0"
attr(fit$metadata, "program")
#> [1] "efron_mix"
substr(fit$metadata$stan_file_sha256, 1, 12)
#> [1] "0ea40bf867ec"
```

If you do not set `seed` in
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
the sampler’s seed is taken from the clock, and it is stored either way.
Of the two Stan programs, which are described below, `efron_mix` is the
default.

The sampler settings are stored as an attribute of `metadata`, and so
are the bayesEfron and cmdstanr versions that made the fit.

``` r

settings <- attr(fit$metadata, "sampler_settings")
str(settings[c("chains", "iter_warmup", "iter_sampling", "adapt_delta",
               "max_treedepth", "init")])
#> List of 6
#>  $ chains       : int 4
#>  $ iter_warmup  : int 1000
#>  $ iter_sampling: int 3000
#>  $ adapt_delta  : num 0.9
#>  $ max_treedepth: int 10
#>  $ init         : num 0.5
attr(fit$metadata, "package_version")
#> [1] "0.4.0"
attr(fit$metadata, "cmdstanr_version")
#> [1] "0.8.0"
```

Another attribute is the grid, stored as
[`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
returned it, together with the method that chose its range. By default
this is `"pilot"`, which fits the model once on a wide range before it
sets the limits of the grid. Because this pilot fit runs with the seed
of the call, the grid is reproduced along with the rest.

``` r

grid <- attr(fit$metadata, "grid")
grid$grid_method
#> [1] "pilot"
grid$support
#> [1] -0.08289621  0.15852607
str(grid$pilot[c("narrowed", "seed", "iter_warmup", "iter_sampling")])
#> List of 4
#>  $ narrowed     : logi FALSE
#>  $ seed         : int 20250807
#>  $ iter_warmup  : int 500
#>  $ iter_sampling: int 500
```

Here the pilot fit used the same seed as the fit and did not narrow the
range.

Finally, `fit$metadata$data_list` contains the data and the spline
basis. If you use the same data, settings, seed and versions on the same
machine, you get the same draws. On another machine or with another
CmdStan version the draws can differ, and the results then agree only to
within Monte Carlo error.

## Stored draws and the size of a fit

By default a fit stores the draws of the spline coefficients, `lambda`
and \\g\\, with five quantities computed from them, and nothing for the
sites. Every site summary is computed from the draws of \\g\\ when a
method needs it, as
[`vignette("model")`](https://joonho112.github.io/bayesEfron/articles/model.md)
describes, so the size of a fit does not grow with the number of sites.

The dimensions of the draws and the `storage` attribute show what the
example fit contains.

``` r

dim(fit$draws)
#> [1] 200   4 114
str(attr(fit$metadata, "storage")[c("site_draws", "iterations", "chains",
                                    "thin")])
#> List of 4
#>  $ site_draws: logi FALSE
#>  $ iterations: int 200
#>  $ chains    : int 4
#>  $ thin      : int 15
```

We made the example fit with `thin = 15`, so it has 800 draws of 114
variables. Without thinning the same fit has 12,000 draws, which take
about 10 megabytes. A fit made with `site_draws = TRUE` also stores five
variables for each site, as version 0.3.0 did, and its 12,000 draws take
about 55 megabytes for the 97 firms and about 697 megabytes for 1,500
sites. These sizes are for the draws alone, and the fitted object as a
whole is a little larger.

Thinning only reduces what is stored. The sampler diagnostics still
describe the whole run, as
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
explains. You can also store the cmdstanr fit object, with
`keep_cmdstan_fit = TRUE`, but it refers to temporary files and is only
useful in the session that created it.

### Fits saved by earlier versions

A fit that was saved by an earlier version of the package can still be
used, although it may lack some of the attributes described above. If
its table of site effects or its `effective_params` was computed under
earlier definitions, the methods compute them again from the stored
draws of \\g\\ and issue a warning, and the saved object is not changed
(see
[`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)).
A saved result of [`summary()`](https://rdrr.io/r/base/summary.html)
contains no draws, so it cannot be updated in this way:
[`print()`](https://rdrr.io/r/base/print.html) throws an error that asks
you to call [`summary()`](https://rdrr.io/r/base/summary.html) again on
the fit. For a fit saved by the first released version, 0.1.0,
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`predict()`](https://rdrr.io/r/stats/predict.html) and `loo()` are not
available and throw an error.

## Two Stan programs

bayesEfron contains two Stan programs. To read them, find the files with
[`system.file()`](https://rdrr.io/r/base/system.file.html):

``` r

system.file("stan", "efron_mix.stan", package = "bayesEfron")
system.file("stan", "efron_re.stan", package = "bayesEfron")
```

The default program is `efron_mix.stan`, and its comments state the
model that
[`vignette("model")`](https://joonho112.github.io/bayesEfron/articles/model.md)
describes. It allows a normal component at each grid point, computes the
likelihood of all the sites as one matrix product, and writes nothing
for the sites.

The other program, `efron_re.stan`, is the program of version 0.3.0. We
have left it as it was, and a test in the package checks that the file
has the same hash as the released one. For a \\g\\ with point masses the
two programs describe the same model, and a test checks that their fits
of the callback data agree within Monte Carlo error.

Each program is compiled the first time it is needed, and the compiled
model is stored in your user cache directory and reused in later
sessions. Its file name includes a hash of the Stan file, the CmdStan
version and the machine architecture, so the program is compiled again
when the Stan file or the CmdStan version changes. See
[`?bayes_efron_compile`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
for where the directory is and how to clear it.

### Reproducing a result of version 0.3.0

To obtain a result of version 0.3.0 again, set `settings = "0.3.0"`. The
fit then runs `efron_re.stan` with the grid method, the 101 grid points
and the site-level draws of that version:

``` r

bayes_efron_fit(raudenbush1985, seed = 1985, settings = "0.3.0")
```

For the same seed this gives the draws of version 0.3.0, provided that
the CmdStan version, the compiler and the machine are the same as well.
Another CmdStan version can change the draws, and also the number of
digits with which they are written.

The current program writes its draws with 12 significant digits, and the
site summaries are computed from the draws of \\g\\ (see above), so
estimates far from zero need no centering. In a fit with
`settings = "0.3.0"`, however, the site-level draws are written with
CmdStan’s default number of digits (8 in recent versions). They lose the
digits in which the sites differ when the estimates are about 10 million
times their standard errors or more, so center such estimates first.

For a fit made with `settings = "0.3.0"`, the summaries follow the
definitions of the current version. To get the interval limits as
version 0.3.0 computed them, use `confint(method = "draws")`.

## The data and the stored results

The vignettes use three data sets, each with a stored fit (`krw_fit`,
`twin_towers_fit` and `raudenbush_fit`).

- `krw_discrimination` contains the callback gaps of 97 firms from the
  experiment of Kline et al. ([2022](#ref-kline2022)). It uses the
  sample of Kline et al. ([2024](#ref-kline2024)), and we recomputed the
  estimates from the replication files of Walters
  ([2024](#ref-walters2024)).
- `twin_towers` is one simulated data set with 1,500 sites from the
  simulation study of Lee and Sui ([2025](#ref-lee2025)), taken from
  that study’s files.
- `raudenbush1985` contains the 19 studies of Raudenbush and Bryk
  ([1985](#ref-raudenbush1985)).

A fourth data set, `star_schools`, has the estimated effects of small
classes in the 79 schools of Tennessee’s Project STAR. The case study
that analyzes it is an article on the package website, and its fit is
not installed with the package.

See each data set’s help page for its source.

We made the fits with `data-raw/example-fits.R` in the source
repository, using the defaults (with 8 spline functions for the
simulated sites) and fixed seeds. The script runs each fit with all
12,000 draws and thins it to 1 draw in 15. It then compares the site
intervals from the thinned draws with the intervals from all the draws,
and it fails with an error if a fit has not converged or an interval
limit moves by more than 2 grid spacings.

For the simulated sites we also compared the probability that a site’s
effect is positive, read on the grid points (`reading = "step"`),
between the thinned draws and all the draws. It differs by at most
0.0013.

No model is fitted when the documentation is built. The vignettes read
these fits and the files in `inst/examples` and `inst/verification`, and
both directories have a README that names the script that made every
file.

## Checks against known answers

The checks in this section are tests in the package’s test suite
(`tests/testthat` in the source repository), so they are run again
whenever the package is checked. The tests that sample need CmdStan.

The tests check the package’s computations in two steps. In the first
step the model’s formulas, written out in plain R, are compared with
results obtained in another way:

- the posterior distribution of a site effect under a normal \\g\\, with
  the conjugate normal result
- the mean and variance of \\g\\, with Gauss-Hermite quadrature
- the posterior mode, with a numerical optimizer

In the second step the quantities that each Stan program writes are
recomputed from the spline coefficients with those R formulas, and the
two are compared.

With a single normal component the default program is the normal model,
so another test samples from it in that case and checks that the result
agrees with the normal model. For the spline basis the package calls
`splines::ns(grid, df = M, intercept = FALSE)`, and a test checks that
it uses the result unchanged. If the columns of the basis and a constant
column are not linearly independent,
[`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
throws an error.

Other tests compare fits with deconvolveR, the implementation by
Narasimhan and Efron ([2020](#ref-narasimhan2020)) of the empirical
Bayes estimator of Efron ([2016](#ref-efron2016)), on three simulated
data sets with equal standard errors and 20, 50 and 100 sites. We do not
expect the two methods to coincide, so the tests only require them to be
close: the posterior means of the site effects within one standard error
of the deconvolveR estimates, with a correlation of at least 0.99, and
the two estimates of \\g\\ for the largest data set within 0.25 in total
variation distance.

These tests show that the pieces are computed correctly. They do not
show that the intervals have the right coverage, which only simulated
data can.

## Interval coverage in simulated data

Whether 90% intervals contain the true effects 90% of the time can only
be checked where the true effects are known. We fitted three sets of
simulated data with the defaults, and the code below collects the stored
results in one table.

``` r

verification <- function(name) {
  read.csv(system.file("verification", name, package = "bayesEfron"),
           check.names = FALSE)
}
panels <- verification("panels-methods.csv")
panels <- panels[panels$method == "bayesEfron", ]
few <- verification("few-studies.csv")
few <- few[few$method == "bayesEfron", ]
equal <- readRDS(system.file("examples", "ebnm_comparison.rds",
                             package = "bayesEfron"))$equal_se$repeated
equal <- equal[equal$method == "bayesEfron" & equal$M == 6, ]

checked <- data.frame(
  effects = c("two groups", "two groups", paste("19 studies,", few$scenario),
              "two groups, equal standard errors"),
  sites = c(panels$K, rep(19, nrow(few)), 1500),
  data_sets = c(panels$panels, few$replications, equal$panels),
  inclusion = c(panels$inclusion, few$inclusion, equal$inclusion),
  se = c(panels$inclusion_se, few$inclusion_se, equal$inclusion_se)
)
checked[c("inclusion", "se")] <- round(checked[c("inclusion", "se")], 3)
checked
#>                             effects sites data_sets inclusion    se
#> 1                        two groups   100       100     0.901 0.004
#> 2                        two groups  1500        20     0.900 0.003
#> 3       19 studies, normal, SD 0.14    19       200     0.904 0.006
#> 4       19 studies, normal, SD 0.31    19       200     0.899 0.005
#> 5              19 studies, fitted g    19       200     0.907 0.006
#> 6 two groups, equal standard errors  1500        20     0.876 0.007
```

In this table, `inclusion` is the proportion of the true effects that
lie inside the 90% intervals of
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
averaged over the data sets. It estimates the coverage of those
intervals, and `se` is its Monte Carlo standard error. There are three
simulation designs.

- In rows 1 and 2, every data set has its own true effects, which we
  simulated from two separate groups: a third of the sites between
  \\-1.7\\ and \\-0.7\\, and two thirds between \\0.7\\ and \\2.7\\.
  These are the populations of the simulation study of Lee and Sui
  ([2025](#ref-lee2025)). The standard errors are unequal, and the
  largest is 9 times the smallest.
- Rows 3 to 5 are the data sets from
  [`vignette("teacher-expectancy")`](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md),
  which have the standard errors of the 19 studies and true effects from
  three distributions.
- In row 6, the true effects of the 20 large data sets in row 2 are each
  estimated with a standard error of 1.

In the five rows with unequal standard errors the coverage is between
0.899 and 0.907, which is within two standard errors of 0.90 in each
row. When every standard error is equal to 1, however, it is 0.876. The
fits in that row converged, and we have not isolated the cause of the
shortfall. See
[`vignette("empirical-bayes")`](https://joonho112.github.io/bayesEfron/articles/empirical-bayes.md)
for what is known about it and what it means in practice.

We also fitted the 20 large data sets in row 2 with 8 spline functions
and with `settings = "0.3.0"`. The next table compares the three
settings.

``` r

convergence <- verification("panels-convergence.csv")
data.frame(setting = convergence$setting, valid = convergence$valid,
           rhat_max = round(convergence$rhat_max, 3),
           ess_bulk_min = round(convergence$ess_bulk_min),
           treedepth_share_max = round(convergence$treedepth_share_max, 3),
           seconds = round(convergence$sampling_seconds_median),
           inclusion = round(convergence$inclusion, 3))
#>              setting valid rhat_max ess_bulk_min treedepth_share_max seconds
#> 1           defaults    20    1.009          679               0.131      36
#> 2              M = 8    20    1.004         1538               0.000      15
#> 3 settings = "0.3.0"     5    1.084           55               0.823     145
#>   inclusion
#> 1     0.900
#> 2     0.893
#> 3     0.900
```

Here `valid` counts the fits that meet three conditions: the largest
R-hat over `alpha`, `lambda` and the mean and standard deviation of
\\g\\ is below 1.01, the smallest bulk effective sample size is at least
400, and there is no divergent transition. The next three columns give
the worst case among the 20 fits, and `seconds` is the median sampling
time.

All 20 fits with the defaults meet the three conditions. Sampling takes
a median of 36 seconds, after 32 seconds for the pilot fit that sets the
range. With `settings = "0.3.0"`, 5 of the 20 fits meet them, and
sampling takes 145 seconds. Sampling is fastest with 8 spline functions,
but the intervals from those fits have a coverage of 0.893, more than
two standard errors below 0.90.
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
uses the same table to discuss what the diagnostics do and do not
measure.

These simulations cover one shape of \\g\\ with many sites and three
shapes with few sites, two patterns of unequal standard errors, and
normal errors with known standard errors. They do not cover numbers of
sites other than 19, 100 and 1,500, standard errors that are themselves
poorly estimated, or effects that are related to their standard errors.
You can find the scripts in `data-raw/benchmark/` in the source
repository, and `inst/verification/README.md` describes the columns of
the tables.

## Relation to the paper

bayesEfron implements the fully Bayesian model of Lee and Sui ([2025,
sec. 4](#ref-lee2025)). The table shows how the paper’s quantities are
written in the package.

| Paper | Package | Meaning |
|:---|:---|:---|
| \\\Theta_i\\ | \\\theta_i\\ in the documentation | True effect of site \\i\\ |
| \\\hat\theta_i\\, \\\sigma_i\\ | `theta_hat`, `sigma` | Estimate and its standard error |
| \\\mathcal{T} = \\\theta_1, \ldots, \theta_m\\\\ | `grid`, with `L` points \\\tau_j\\ | Grid on which \\g\\ is defined |
| \\Q\\, an \\m \times p\\ matrix | `B`, an `L` by `M` matrix | Spline basis on the grid |
| \\\alpha\\, \\\lambda\\ | `alpha`, `lambda` | Spline coefficients and their prior precision |
| \\g_j\\ | `exp(log_g)` | Probability of grid point \\j\\ |
| \\\pi\_{ij}\\ | \\w\_{ij}\\ in [`vignette("model")`](https://joonho112.github.io/bayesEfron/articles/model.md); not stored | Posterior probability of grid point \\j\\ for site \\i\\, given \\\alpha\\ |

In three places the package’s defaults differ from the paper, and each
can be changed back to the paper’s choice.

1.  In the code of the paper’s application, the grid covers the range of
    the estimates extended by half its width on each side, which is
    `grid_method = "paper_realdata"`. The default method, `"pilot"`,
    starts from that range and narrows it where the pilot fit shows that
    \\g\\ occupies only a part of it, but it never widens it. See
    [`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)
    for details.
2.  The number of grid points is 101 unless the smallest standard error
    or the spread of the estimates is small compared with the grid
    spacing, in which case more points are used. Set `L = 101` to keep
    101 points.
3.  The paper obtains the site intervals as quantiles of posterior draws
    of the effects (its Appendix B.4), and so did version 0.3.0 of the
    package. We now compute them from the posterior distribution itself,
    with the probability of each grid point spread over its grid cell.
    In
    [`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
    `reading = "step"` leaves the probabilities on the points, and a fit
    made with `site_draws = TRUE` gives the intervals from draws with
    `confint(method = "draws")`.

To restore the first two at once, use `settings = "0.3.0"`, which also
stores the site-level draws. The other two grid methods,
`"paper_simulation"` and `"paper_sensitivity"`, are the methods that the
code accompanying the paper,
<https://github.com/joonho112/fully-bayes-efron-prior>, uses for its
simulation study and its sensitivity analysis. Where the paper’s text
and its code differ in a detail, the package follows the code. For
example, the paper’s Stan appendix shows a grid that extends 0.5 units
beyond the estimates on each side, whereas the simulation code builds
the grid from the true effects.

[`confint()`](https://rdrr.io/r/stats/confint.html) can also return the
shortest interval with posterior probability `level`
(`interval = "hpd"`). For point masses it is the shortest interval after
the probability of each grid point has been spread evenly over one grid
spacing around it, and where a range of intervals is equally short, the
middle one is returned. For normal components the limits are moved until
no nearby interval with that probability is shorter, and the length of
the result is within one hundredth of the length of the shortest
interval. In both cases the interval has probability `level` exactly.
With `method = "draws"`, and always for `type = "g"`, it is the shortest
interval that contains at least that proportion of the draws.

The package also has features that are not in the paper:

- the normal components
- the `"pilot"` method for the range, and what
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  reports about the grid
- the effective number of parameters
- the log likelihood of each site and the cross-validation built on it
- [`predict()`](https://rdrr.io/r/stats/predict.html)
- the plots

The table of site effects follows the paper’s definitions ([Lee and Sui
2025, sec. 4.3](#ref-lee2025) and Appendix A): the posterior mean, the
posterior standard deviation with both parts of the variance, the
posterior mode on the grid, and an interval of the stated probability,
which the package takes with equal tails.

To cite the package, please cite the paper. You can get the reference
with `citation("bayesEfron")`.

## References

Efron, Bradley. 2016. “Empirical Bayes Deconvolution Estimates.”
*Biometrika* 103 (1): 1–20. <https://doi.org/10.1093/biomet/asv068>.

Kline, Patrick, Evan K. Rose, and Christopher R. Walters. 2022.
“Systemic Discrimination Among Large U.S. Employers.” *The Quarterly
Journal of Economics* 137 (4): 1963–2036.
<https://doi.org/10.1093/qje/qjac024>.

Kline, Patrick, Evan K. Rose, and Christopher R. Walters. 2024. “A
Discrimination Report Card.” *American Economic Review* 114 (8):
2472–525. <https://doi.org/10.1257/aer.20230700>.

Lee, JoonHo, and Daihe Sui. 2025. “Fully Bayesian Inference for
Meta-Analytic Deconvolution Using Efron’s Log-Spline Prior.”
*Mathematics* 13 (16): 2639. <https://doi.org/10.3390/math13162639>.

Narasimhan, Balasubramanian, and Bradley Efron. 2020. “deconvolveR: A
G-Modeling Program for Deconvolution and Empirical Bayes Estimation.”
*Journal of Statistical Software* 94 (11): 1–20.
<https://doi.org/10.18637/jss.v094.i11>.

Raudenbush, Stephen W., and Anthony S. Bryk. 1985. “Empirical Bayes
Meta-Analysis.” *Journal of Educational Statistics* 10 (2): 75–98.
<https://doi.org/10.3102/10769986010002075>.

Walters, Christopher R. 2024. “Empirical Bayes Methods in Labor
Economics.” In *Handbook of Labor Economics*, vol. 5. Elsevier.
<https://doi.org/10.1016/bs.heslab.2024.11.001>.
