# Leave-one-out cross-validation and WAIC for a fitted model

`loo()` estimates how well the fitted model predicts the estimate of a
site that is left out. It applies Pareto smoothed importance sampling
(Vehtari, Gelman and Gabry, 2017) to the log likelihood of each site,
which is computed from the draws of \\g\\. Use the result to compare
fits of the same data that differ in the grid or in the number of spline
functions;
[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
makes that comparison in one call.

`waic()` computes the widely applicable information criterion from the
same log likelihood. We recommend `loo()`, because its Pareto k values
show for each site whether the approximation can be relied on, whereas
`waic()` issues only a general warning when some of its terms are large.

## Usage

``` r
# S3 method for class 'bef_fit'
loo(x, ...)

# S3 method for class 'bef_fit'
waic(x, ...)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- ...:

  Passed to [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)
  or [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html).

## Value

An object of class `psis_loo` from
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html), or of class
`waic` from
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html).

## Details

Both methods need the loo package. The relative efficiency that
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) uses for its
Monte Carlo standard errors is computed from the draws with
[`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html).

## References

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413-1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

## See also

[`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md),
[`plot_loo()`](https://joonho112.github.io/bayesEfron/reference/plot_loo.md)

## Examples

``` r
fit_loo <- loo::loo(krw_fit)
fit_loo
#> 
#> Computed from 800 by 97 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo    227.4  7.7
#> p_loo         2.3  0.7
#> looic      -454.7 15.4
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [1.0, 1.1]).
#> 
#> All Pareto k estimates are good (k < 0.66).
#> See help('pareto-k-diagnostic') for details.
range(fit_loo$diagnostics$pareto_k)
#> [1] -0.1887404  0.3386994
```
