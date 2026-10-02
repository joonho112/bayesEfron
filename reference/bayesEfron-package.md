# bayesEfron: Fully Bayesian Inference for the Empirical-Bayes Deconvolution Problem

Estimates the distribution of true effects from a set of estimates with
known standard errors, as in a meta-analysis or a multisite trial, using
the fully Bayesian version of Efron's log-spline prior described in Lee
and Sui (2025)
[doi:10.3390/math13162639](https://doi.org/10.3390/math13162639) . The
fit returns posterior draws and summaries for each effect and for the
effect distribution, with sampler diagnostics, plots, and model
comparison by approximate leave-one-out cross-validation. The model is
written in 'Stan' and is fitted with 'CmdStan' through the 'cmdstanr'
package.

## Details

For sites \\i = 1, \ldots, K\\ the model is \$\$\hat\theta_i \mid
\theta_i \sim N(\theta_i, \sigma_i^2), \qquad \theta_i \sim g,\$\$ where
the standard errors \\\sigma_i\\ are taken as known and may differ
between sites. The distribution \\g\\ is represented on a grid, with log
probabilities that are a combination of natural cubic spline functions;
the spline coefficients have a normal prior whose precision has a
half-Cauchy prior.
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
gives the details.

## Main functions

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  takes the estimates and standard errors out of a list, a data frame or
  a
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  object and checks them.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  fits the model and returns a
  [bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
  object, which has [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and other methods.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  returns the convergence diagnostics of a fit.

- [plot()](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  draws the site effects, the estimated distribution of effects, the
  diagnostics or the Pareto k values.

- [predict()](https://joonho112.github.io/bayesEfron/reference/predict.bef_fit_re.md)
  draws the effect or the estimate of a new site.

- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  builds the grid and the spline basis, and
  [`compare_efron_grids()`](https://joonho112.github.io/bayesEfron/reference/compare_efron_grids.md)
  compares several of them by leave-one-out cross-validation.

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  compiles the Stan model before the first fit, and
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  removes the compiled model.

The examples use the data set
[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
and the fit
[raudenbush_fit](https://joonho112.github.io/bayesEfron/reference/raudenbush_fit.md).
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)
is the place to start.

## Funding

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama. The opinions expressed are those of the authors and do not
represent views of the Institute or the U.S. Department of Education.

## See also

Useful links:

- <https://github.com/joonho112/bayesEfron>

- <https://joonho112.github.io/bayesEfron/>

- Report bugs at <https://github.com/joonho112/bayesEfron/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
