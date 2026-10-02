# Fit the Efron log-spline model to estimates and standard errors

A meta-analysis or a multisite trial gives an estimate \\\hat\theta_i\\
and a standard error \\\sigma_i\\ for each of \\K\\ sites.
`bayes_efron_fit()` treats the true effects \\\theta_i\\ as draws from
an unknown distribution \\g\\ and estimates \\g\\ and each \\\theta_i\\
together. Unlike the usual random-effects model, it does not take \\g\\
to be normal: \\g\\ can be skewed or have more than one mode, and how
far each estimate is shrunk depends on that shape.

The function returns posterior draws of \\g\\. For each site and each
draw of \\g\\ it also returns the mean, standard deviation and mode of
\\\theta_i\\ given that \\g\\, and one draw of \\\theta_i\\. The model
is fitted by Markov chain Monte Carlo with CmdStan, so the cmdstanr
package and CmdStan must be installed. The first call compiles the Stan
model, which can take up to a minute; later calls reuse the compiled
model (see
[`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)).

## Usage

``` r
bayes_efron_fit(
  theta_hat,
  sigma,
  ...,
  grid_method = c("paper_realdata", "paper_simulation", "paper_sensitivity",
    "kl_target_experimental"),
  L = 101L,
  expansion = 0.5,
  M = 6L,
  theta_true = NULL,
  bound_expansion = NULL,
  model_family = "RE",
  chains = 4L,
  parallel_chains = chains,
  iter_warmup = 1000L,
  iter_sampling = 3000L,
  adapt_delta = 0.9,
  max_treedepth = 10L,
  seed = NULL,
  keep_cmdstan_fit = FALSE,
  store_grid_quantities = FALSE
)
```

## Arguments

- theta_hat:

  Numeric vector of effect estimates, one for each site, all on the same
  scale (standardized mean differences or log odds ratios, for example).
  At least five are needed.

- sigma:

  Numeric vector of the standard errors of the estimates, positive and
  of the same length as `theta_hat`.

- ...:

  These dots are for future extensions and must be empty.

- grid_method:

  Character string, the rule for the end points of the grid:
  `"paper_realdata"` (the default), `"paper_simulation"`,
  `"paper_sensitivity"` or `"kl_target_experimental"`. The Details of
  [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  describe them.

- L:

  Whole number between 51 and 300, the number of grid points. Defaults
  to `101L`.

- expansion:

  Number between 0 and 5. With `"paper_realdata"` and
  `"kl_target_experimental"`, the range of `theta_hat` is widened on
  each side by `expansion` times its width. Defaults to `0.5`.

- M:

  Whole number between 3 and 10, the number of spline functions.
  Defaults to `6L`.

- theta_true:

  Numeric vector of the true effects, of the same length as `theta_hat`.
  Needed by `"paper_simulation"` and `"paper_sensitivity"`, which are
  for simulations.

- bound_expansion:

  Number greater than 0 and at most 5. With `"paper_sensitivity"`, the
  range of `theta_true` is widened on each side by `bound_expansion`
  times its width. `NULL` (the default) means 0.5.

- model_family:

  Character string. `"RE"`, the model described in Details, is the only
  choice.

- chains:

  Whole number between 1 and 16, the number of Markov chains. Defaults
  to `4L`.

- parallel_chains:

  Whole number between 1 and `chains`, the number of chains to run at
  the same time. Defaults to `chains`.

- iter_warmup:

  Whole number, the number of warmup iterations in each chain. Defaults
  to `1000L`.

- iter_sampling:

  Whole number of at least 1, the number of iterations in each chain
  that are kept. Defaults to `3000L`.

- adapt_delta:

  Number between 0 and 1, the target acceptance rate of the sampler.
  Defaults to `0.9`. A value nearer 1 gives smaller steps, which can
  remove divergent transitions at the cost of a longer run.

- max_treedepth:

  Whole number between 1 and 20, the maximum tree depth of the sampler.
  Defaults to `10L`. Raise it if
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  reports many iterations that reached the maximum.

- seed:

  Whole number, the seed of the sampler, or `NULL` (the default), in
  which case a seed is taken from the clock. The seed that was used is
  stored in the result.

- keep_cmdstan_fit:

  Logical. If `TRUE`, the
  [`cmdstanr::CmdStanMCMC`](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)
  object is kept as the element `cmdstan_fit` of the result. It refers
  to temporary files, so it is of use only in the session that made it.
  Defaults to `FALSE`.

- store_grid_quantities:

  Logical. If `TRUE`, the draws also hold the grid probabilities `g` and
  their unnormalized logarithms `log_w`, which adds `2 * L` variables.
  Defaults to `FALSE`; the methods compute both from the other draws
  when they need them.

## Value

An object of class `bef_fit_re`, which inherits from `bef_fit`. Its
elements are `draws`, the posterior draws as a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html);
`posterior`, the draws of the mean, variance and standard deviation of
\\g\\ and of two summaries of the fit; and `metadata`, the data, the
settings and summaries of the draws.
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
describes the elements and the methods.

## Details

For sites \\i = 1, \ldots, K\\, \$\$\hat\theta_i \mid \theta_i \sim
N(\theta_i, \sigma_i^2), \qquad \theta_i \sim g.\$\$ The standard errors
\\\sigma_i\\ are taken as known and may differ between sites. The
distribution \\g\\ has probabilities \\g_1, \ldots, g_L\\ on a grid of
points \\\tau_1 \< \cdots \< \tau_L\\, with \$\$\log g_j = B_j^\top
\alpha - \log \sum\_{l=1}^{L} \exp(B_l^\top \alpha),\$\$ where \\B_j\\
holds the values of the \\M\\ natural cubic spline functions at
\\\tau_j\\. This is the log-spline family of Efron (2016). The
coefficients have independent priors \\\alpha_m \mid \lambda \sim N(0,
1/\lambda)\\, which draw \\g\\ toward the uniform distribution on the
grid, and the precision has the prior \\\lambda \sim
\mbox{half-Cauchy}(0, 5)\\, so that the strength of that pull is
estimated along with the coefficients. Lee and Sui (2025) describe the
model and compare it with the empirical Bayes estimator, in which
\\\alpha\\ is fixed at an estimate.

The grid and the spline basis are those of
[`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md),
and `grid_method`, `L`, `expansion` and `M` control them. The defaults,
101 grid points and six spline functions on a grid that extends beyond
the estimates by half their range on each side, follow the application
in Lee and Sui (2025, Section 7.2) and its replication code.
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)
discusses when to change them.

Only \\\alpha\\ and \\\lambda\\ are sampled. Given \\\alpha\\, the
posterior distribution of each \\\theta_i\\ on the grid has a closed
form, and its mean, standard deviation and mode and one draw from it are
computed for every posterior draw of \\\alpha\\.

Each chain starts from spline coefficients drawn uniformly between -0.5
and 0.5 and from a precision whose logarithm is drawn in the same way, a
narrower interval than CmdStan's default of -2 to 2. Pass a `seed` to
make a fit reproducible: with the same data, settings, versions of the
packages and of CmdStan, and machine, the draws are the same. The seed
and the settings of the sampler are stored with the result;
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
says where.

The function gives a warning if there are divergent transitions and, for
a fit with at least 400 draws in all, if an R-hat value exceeds 1.05;
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
returns all the diagnostics.

CmdStan writes the draws with a limited number of significant digits,
eight in recent versions. That is ample for effect sizes. Estimates that
are very large compared with their standard errors, by a factor of
100,000 or more, should be centered before they are passed, because the
saved draws would otherwise lose the digits in which the sites differ.

## References

Efron, B. (2016). Empirical Bayes deconvolution estimates. *Biometrika*,
103(1), 1-20.
[doi:10.1093/biomet/asv068](https://doi.org/10.1093/biomet/asv068)

Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
2639. [doi:10.3390/math13162639](https://doi.org/10.3390/math13162639)

## See also

[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md),
[`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md),
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md),
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
[`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md),
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)

## Examples

``` r
# The call that produced `raudenbush_fit`. It needs CmdStan.
if (FALSE) { # \dontrun{
fit <- bayes_efron_fit(
  theta_hat = raudenbush1985$yi,
  sigma = raudenbush1985$sei,
  seed = 1985
)
} # }

# The stored result of that call.
fit <- raudenbush_fit
fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 19
#> Grid method: paper_realdata
#> Runtime: 1.351 sec
#> Diagnostics: Rhat 1.007; ESS bulk 676.6; ESS tail 888.8; divergences 0; treedepth hits 1
#> Warnings: max_treedepth
#> Use summary() for posterior summaries.
summary(fit)
#> <summary.bef_fit>
#> 
#> Prior g:
#>   mean: 0.1145
#>   var:  0.1094
#>   sd:   0.3009
#> 
#> Diagnostics:
#>   Rhat (max):        1.007
#>   ESS bulk (min):    676.6
#>   ESS tail (min):    888.8
#>   Divergences:       0
#>   Treedepth hits:    1
#>   Effective params:  mean 10.29, sd 2.84
#>   Log marginal lik.: mean -4.047, sd 1.439
#>   Runtime:           1.351 sec
#>   Warnings:          max_treedepth
#> 
#> Theta summary:
#>   site      mean        sd     lower     upper       map
#>   1      0.03805   0.09832     -0.11      0.19      0.04
#>   2       0.0874    0.1121     -0.08      0.28      0.07
#>   3     -0.03442    0.1188     -0.23      0.16     -0.02
#>   4        0.635    0.4979      0.04      1.48      0.19
#>   5       0.1063     0.186     -0.14    0.4015      0.07
#>   6     -0.02282   0.08604     -0.17       0.1     -0.02
#>   14      0.1073    0.1658     -0.11       0.4      0.07
#>   15    -0.05756    0.1177     -0.29    0.1015     -0.05
#>   16    0.001513    0.1169      -0.2      0.19      0.01
#>   17      0.1979    0.1221      0.01       0.4      0.19
#>   18       0.063   0.08101     -0.08      0.19      0.07
#>   19   -0.0006152    0.1197      -0.2      0.19      0.01
#>   ... 7 sites omitted ...
confint(fit, parm = c(4, 10))
#>   site lower upper     point
#> 1    4  0.04  1.48 0.6349621
#> 2   10  0.07  1.03 0.4652042
plot(fit, type = "g")

```
