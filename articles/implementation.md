# Implementation and checks

Someone who reports results from this package may be asked whether the
fit can be reproduced and what evidence there is that the intervals mean
what they say.

## Reproducing a fit

A fitted object records what produced its draws.

``` r

fit$metadata$seed
#> [1] 1985
fit$metadata$cmdstan_version
#> [1] "2.38.0"
substr(fit$metadata$stan_file_sha256, 1, 12)
#> [1] "95fa186dda55"
```

`seed` is the seed of the sampler. If none is given to
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
one is taken from the clock, and it is stored either way. The last value
is the beginning of a hash of the Stan file. The settings of the sampler
and the versions of bayesEfron and cmdstanr are attributes of
`metadata`.

``` r

str(attr(fit$metadata, "sampler_settings"))
#> List of 9
#>  $ chains         : int 4
#>  $ parallel_chains: int 4
#>  $ iter_warmup    : int 1000
#>  $ iter_sampling  : int 3000
#>  $ adapt_delta    : num 0.9
#>  $ max_treedepth  : int 10
#>  $ seed           : int 1985
#>  $ refresh        : int 0
#>  $ init           : num 0.5
attr(fit$metadata, "package_version")
#> [1] "0.3.0"
attr(fit$metadata, "cmdstanr_version")
#> [1] "0.8.0"
```

The data, the grid and the spline basis are in `fit$metadata$data_list`.
With the same data, settings, seed and versions on the same machine, the
draws are the same. Between machines or versions of CmdStan they can
differ, and the results then agree to within Monte Carlo error only.

## The Stan model

The model is one Stan file, which can be read at

``` r

system.file("stan", "efron_re.stan", package = "bayesEfron")
```

Its comments state the model, which
[`vignette("model")`](https://joonho112.github.io/bayesEfron/articles/model.md)
describes. The file is compiled the first time a model is fitted. The
compiled model is kept in the user’s cache directory and used again in
later sessions, and it is recompiled when the Stan file or the version
of CmdStan changes.
[`?bayes_efron_compile`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
says where the directory is and how to clear it.

## Checks against known answers

The tests of the package check its computations in two steps. First the
formulas of the model are written out in plain R and compared with
results obtained in another way: the posterior distribution of a site
effect under a normal \\g\\ with the conjugate normal result, and the
mean and variance of \\g\\ with Gauss-Hermite quadrature. Then the
quantities that the Stan model wrote for the stored fit `raudenbush_fit`
are recomputed with those R formulas from the draws of the spline
coefficients and compared: the probabilities of \\g\\, its mean and
variance, the posterior mean, standard deviation and mode of each site,
the log likelihood and the effective number of parameters.

Other tests compare fits with deconvolveR, the implementation by
Narasimhan and Efron ([2020](#ref-narasimhan2020)) of the empirical
Bayes estimator of Efron ([2016](#ref-efron2016)), on three simulated
data sets with equal standard errors and 20, 50 and 100 sites. The two
methods are not expected to coincide. The tests require that the
posterior means of the site effects lie within one standard error of the
estimates of deconvolveR, with a correlation of at least 0.99, and that
the two estimates of \\g\\ for the largest data set differ by at most
0.25 in total variation distance.

These tests show that the pieces are computed correctly. They do not
show that the intervals have the right coverage.

## A benchmark with known effects

The benchmark repeats, with the package and with Markov chain Monte
Carlo, a check that Lee and Sui ([2025](#ref-lee2025), Appendix E) made
with a variational approximation. The true effects are the 1,500 values
of their simulation study, in two separate groups: 500 between \\-1.7\\
and \\-0.7\\ and 1,000 between \\0.7\\ and \\2.7\\. The standard errors
are unequal, the largest being 9 times the smallest. Data sets of 50,
100, 200, 500 and 1,500 sites are fitted, 20 of each size, with the
default sampler settings, the default rule for the ends of the grid and
the numbers of grid points used in the paper (51, 71, 81, 101 and 101).

The 20 data sets of 50 to 500 sites are stratified random samples of the
sites of one simulated panel of 1,500. They differ in the sites that
they hold, and they overlap, more so the larger they are. The 20 data
sets of 1,500 sites hold the same true effects, each with its own
assignment of the standard errors to the sites and its own errors. The
coverage of a fit is the share of sites whose true effect lies inside
the 90% interval that
[`confint()`](https://rdrr.io/r/stats/confint.html) returns.

``` r

benchmark <- utils::read.csv(
  system.file("verification", "benchmark-coverage.csv", package = "bayesEfron"),
  check.names = FALSE
)
shown <- benchmark[, c("K", "n_replications", "coverage", "coverage_se",
                       "max_rhat", "n_rhat_above_1.05", "min_ess_bulk",
                       "total_divergences", "total_treedepth_hits")]
shown[c("coverage", "coverage_se", "max_rhat")] <-
  round(shown[c("coverage", "coverage_se", "max_rhat")], 3)
shown$min_ess_bulk <- round(shown$min_ess_bulk)
shown
#>      K n_replications coverage coverage_se max_rhat n_rhat_above_1.05
#> 1   50             20    0.891       0.010    1.019                 0
#> 2  100             20    0.891       0.007    1.020                 0
#> 3  200             20    0.891       0.006    1.008                 0
#> 4  500             20    0.891       0.004    1.008                 0
#> 5 1500             20    0.876       0.004    2.352                 3
#>   min_ess_bulk total_divergences total_treedepth_hits
#> 1          362                 0                    5
#> 2          229                 0                 2235
#> 3          400                 0                 9147
#> 4         1006                 0                36411
#> 5            5                 0               167442
```

The table was computed on 2026-10-01 with CmdStan 2.38.0. `coverage_se`
is the standard deviation of the coverage between data sets divided by
the square root of their number. It measures the variation between the
20 data sets of one size and nothing more: the data sets of up to 500
sites are samples of the sites of one simulated panel, and those of
1,500 sites share one set of true effects, so the variation that another
panel or another set of true effects would bring is not in it. The last
five columns describe convergence: the largest R-hat and the smallest
bulk effective sample size over the fits, the number of fits with an
R-hat above 1.05, and the total numbers of divergent transitions and of
iterations that reached the maximum tree depth.

With up to 500 sites the coverage is about 0.89 at every size, a little
below the nominal 0.90. For the same sizes Lee and Sui
([2025](#ref-lee2025)) report 88.7% to 89.6% with the variational
approximation. None of the 80 fits has an R-hat above 1.05, although at
50 and 100 sites the largest R-hat exceeds 1.01 or the smallest bulk
effective sample size is below 400, so some of those fits would have
drawn a warning from
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md).

For 1,500 sites the coverage is 0.876, but 3 of the 20 fits did not
converge with the default settings: the largest R-hat is 2.35, the
smallest bulk effective sample size is 5, and 167,442 of the 240,000
iterations reached the maximum tree depth. The row therefore describes
what the default settings give at that size and is not evidence about
the coverage of a converged fit.

We refitted the data sets of the fits that did not converge in two ways,
with a narrower grid and with a larger maximum tree depth.

``` r

remedies <- utils::read.csv(
  system.file("verification", "benchmark-remedies.csv", package = "bayesEfron")
)
shown <- remedies[, c("replication", "setting", "max_rhat", "min_ess_bulk",
                      "treedepth_hits", "coverage")]
shown$max_rhat <- round(shown$max_rhat, 2)
shown$min_ess_bulk <- round(shown$min_ess_bulk)
shown$coverage <- round(shown$coverage, 3)
shown
#>   replication            setting max_rhat min_ess_bulk treedepth_hits coverage
#> 1          10            default     1.14           20          10166    0.896
#> 2          10   expansion = 0.05     1.00         1628           3088    0.897
#> 3          10 max_treedepth = 12     1.00          829           5554    0.897
#> 4          13            default     2.35            5          10509    0.859
#> 5          13   expansion = 0.05     1.04          149           9849    0.878
#> 6          13 max_treedepth = 12     1.25           12          10055    0.861
#> 7          16            default     1.18           19          10393    0.845
#> 8          16   expansion = 0.05     1.02          201           7989    0.857
#> 9          16 max_treedepth = 12     1.02          122           9329    0.839
```

With a grid that ends just beyond the estimates (`expansion = 0.05` in
place of 0.5), R-hat was at most 1.04 and no fit failed the check,
although 2 of the 3 still had an effective sample size below 400. With
the default grid and a maximum tree depth of 12, 1 of the 3 fits still
had an R-hat above 1.05, and each took several times as long.

In these data many of the standard errors are large, so the estimates
spread over a much wider range than the true effects, and the default
grid, which extends beyond the estimates by half their range on each
side, is wider than it needs to be. This cannot be made a general rule,
however. When the standard errors are small, the estimates spread little
more than the true effects, and room beyond them is needed
([`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)).

The benchmark covers one shape of \\g\\ and one pattern of standard
errors. It says nothing about data whose standard errors are much larger
relative to the spread of the effects, or about fewer than 50 sites. The
scripts and the data are in `data-raw/benchmark/` in the source
repository.

## Relation to the paper

The package implements the fully Bayesian model of Lee and Sui ([2025,
sec. 4](#ref-lee2025)). It writes the quantities of the paper as
follows.

| Paper | Package | Meaning |
|:---|:---|:---|
| \\\Theta_i\\ | \\\theta_i\\ in the documentation | True effect of site \\i\\ |
| \\\hat\theta_i\\, \\\sigma_i\\ | `theta_hat`, `sigma` | Estimate and its standard error |
| \\\mathcal{T} = \\\theta_1, \ldots, \theta_m\\\\ | `grid`, with `L` points \\\tau_j\\ | Grid on which \\g\\ is defined |
| \\Q\\, an \\m \times p\\ matrix | `B`, an `L` by `M` matrix | Spline basis on the grid |
| \\\alpha\\, \\\lambda\\ | `alpha`, `lambda` | Spline coefficients and their prior precision |
| \\g_j\\ | `exp(log_g)` | Probability of grid point \\j\\ |
| \\\pi\_{ij}\\ | \\w\_{ij}\\ in [`vignette("model")`](https://joonho112.github.io/bayesEfron/articles/model.md); not stored | Posterior probability of grid point \\j\\ for site \\i\\, given \\\alpha\\ |

The grid methods `"paper_realdata"`, `"paper_simulation"` and
`"paper_sensitivity"` are the rules of the code that accompanies the
paper, <https://github.com/joonho112/fully-bayes-efron-prior>, for the
application, the simulation study and the sensitivity analysis. Where
the text of the paper and its code differ in a detail, the package
follows the code. The Stan appendix of the paper, for example, shows a
grid that extends 0.5 beyond the estimates, whereas the simulation code
builds it from the true effects.

Not in the paper are the effective number of parameters, the log
likelihood of each site with the cross-validation built on it,
[`predict()`](https://rdrr.io/r/stats/predict.html) and the experimental
grid method. The table of site effects follows the definitions of the
paper ([Lee and Sui 2025, sec. 4.3](#ref-lee2025) and Appendix A): the
posterior mean, the posterior standard deviation with both parts of the
variance, the posterior mode on the grid and an equal-tailed interval.

To cite the package, cite the paper; `citation("bayesEfron")` gives the
reference.

## References

Efron, Bradley. 2016. “Empirical Bayes Deconvolution Estimates.”
*Biometrika* 103 (1): 1–20. <https://doi.org/10.1093/biomet/asv068>.

Lee, JoonHo, and Daihe Sui. 2025. “Fully Bayesian Inference for
Meta-Analytic Deconvolution Using Efron’s Log-Spline Prior.”
*Mathematics* 13 (16): 2639. <https://doi.org/10.3390/math13162639>.

Narasimhan, Balasubramanian, and Bradley Efron. 2020. “deconvolveR: A
G-Modeling Program for Deconvolution and Empirical Bayes Estimation.”
*Journal of Statistical Software* 94 (11): 1–20.
<https://doi.org/10.18637/jss.v094.i11>.
