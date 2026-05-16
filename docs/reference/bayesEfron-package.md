# bayesEfron: Fully Bayesian Inference for the Empirical-Bayes Deconvolution Problem

`bayesEfron` fits the fully Bayesian Efron log-spline prior for
univariate random-effects meta-analytic deconvolution with
heteroscedastic within-study standard errors. The package is for the
applied meta-analyst who has a vector of study-level effect estimates
and their standard errors and wants two related quantities: posterior
summaries of the latent site effects, and a continuous, smoothly
regularised estimate of the underlying mixing distribution that
generated those effects. Both objects are returned from a single fit,
with calibrated credible intervals derived from the same posterior
draws.

## Details

For sites \\i = 1, \ldots, K\\, \$\$\hat\theta_i \mid \theta_i \sim
\mathcal{N}(\theta_i, \sigma_i^2), \qquad \theta_i \mid g \sim g,\$\$
where \\g\\ is the unknown mixing distribution. The package represents
\\g\\ as a discrete distribution on a fixed grid of length \\L\\
(default \\L = 101\\) whose log-density is a linear combination of \\M\\
natural-cubic-spline basis functions (default \\M = 6\\). The spline
coefficients receive a weakly informative Gaussian prior whose precision
\\\lambda\\ is in turn assigned a half-Cauchy hyperprior, so the
smoothness of the deconvolved density is itself part of what the sampler
learns.

Posterior computation is delegated to CmdStan via the `cmdstanr`
package, with two-tier compile caching so that repeated fits with the
same model code reuse the compiled binary. The within-study standard
errors enter exactly: each site's likelihood contribution uses its own
\\\sigma_i\\, with no homoscedasticity assumption.

## User-facing entry points

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md):

  Primary fitting entry. Takes per-site effect estimates and standard
  errors, runs the eight-stage pipeline (validation, grid construction,
  Stan-data preparation, cache-backed model retrieval, sampling, draw
  extraction, postprocessing, assembly), and returns a fitted
  `bef_fit_re` object carrying posterior draws, the deconvolved density,
  and posterior summaries of the site effects.

- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md):

  Constructs the discrete support of \\g\\ and the natural-cubic-spline
  basis evaluated on it. Four recipes: paper real-data, paper
  simulation, paper sensitivity, and an experimental KL-target recipe.

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md):

  Input adapter. Converts a plain list of `theta_hat`/`sigma`, an
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  object, or an existing `bef_data` object into the canonical input
  class.

- [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md):

  Diagnostic producer. Returns a `bef_diagnostic` object summarising
  R-hat, effective sample sizes, divergent transitions, max-treedepth
  saturations, and other sampler-health quantities.

## Where to start

New users should begin with the applied vignette *A1 · Getting started*,
which walks through a minimal end-to-end fit on a five-site toy dataset.
The companion methodological vignette *M1 · The empirical-Bayes
deconvolution problem* develops the statistical background, motivates
the log-spline prior, and explains how the deconvolved density \\\hat
g\\ relates to the posterior site effects \\\theta_i\\.

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

**Maintainer**: JoonHo Lee <joonho112@users.noreply.github.com>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
