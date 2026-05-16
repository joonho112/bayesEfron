# Changelog

## bayesEfron 0.1.0

Initial public snapshot of the random-effects core for fully Bayesian
Efron log-spline deconvolution in univariate meta-analysis with
heteroscedastic standard errors.

### User-facing functionality

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  is the main fitting entry point: it takes per-site effect estimates
  and within-study standard errors, samples the log-spline mixing
  distribution $`g(\theta)`$ jointly with the per-site latent effects
  $`\theta_i`$, and returns a `bef_fit_re` object carrying posterior
  draws, summaries, generated quantities, sampler diagnostics, and
  reproducibility metadata.
- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  constructs the discrete support and natural cubic-spline basis for
  $`g(\theta)`$ and ships four named recipes — `paper_realdata`,
  `paper_simulation`, `paper_sensitivity`, and the experimental
  `kl_target_experimental` — each recording its formula and source
  attribution on the returned object.
- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  is the input converter; it accepts existing `bef_data` objects, plain
  lists of `(theta_hat, sigma)`, and
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  output, so a `metafor` aggregate-data workflow flows directly into a
  bayesEfron fit.
- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  pre-warms the CmdStan model cache and runs a short post-compile smoke
  check, and
  [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_clear_cache.md)
  provides idempotent cleanup at four scopes (`lock_only`, `session`,
  `compiled_models`, `all`) for managing the on-disk cache between
  sessions or across CmdStan upgrades.
- Fitted random-effects objects support the standard model-object
  surface: [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html), and
  [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html),
  plus
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
  for a tidy `bef_diagnostic` view of sampler convergence and posterior
  predictive flags.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for fitted
  random-effects objects offers caterpillar, prior-summary, single-fit
  sensitivity, and diagnostic views, with both a base-graphics default
  and an optional `ggplot2` backend (`BAYESEFRON_NO_GGPLOT2=1` forces
  base).
- [`format()`](https://rdrr.io/r/base/format.html) and
  [`print()`](https://rdrr.io/r/base/print.html) methods for `bef_fit`,
  `summary.bef_fit`, `bef_data`, and `bef_diagnostic` use `cli`
  formatting when available and fall back to plain base output under
  `BAYESEFRON_NO_CLI=1`, so the same code prints cleanly in both
  interactive and log-capture contexts.

### Release status and verification

- v0.1 scope is intentionally restricted to the random-effects model:
  non-`RE` model families, grouped designs, and within-cluster
  correlations are rejected at the input boundary rather than silently
  reinterpreted.
- The v0.1.0 release gate is the Tier 3 calibration of `theta_rep`
  interval coverage on the Lee–Sui benchmark, with all twenty
  replications run per K and the aggregate coverage required to fall
  inside the pre-specified acceptance interval `[0.87, 0.92]` at K
  $`\in`$ {50, 100, 200, 500, 1500}. The K = 1500 result was the thin
  worst case, with aggregate coverage of 0.876433333333333 — inside
  `[0.87, 0.92]` and recorded as the documented release evidence.
- Tier 3 live verification is environment-variable gated.
  `BAYESEFRON_RUN_LIVE=1` opts into the smoke path,
  `BAYESEFRON_RUN_FULL_LIVE=1` opts into the full five-K
  twenty-replication release path, and a fresh full-live refit without
  the accepted replay matrix additionally requires
  `BAYESEFRON_TIER3_OK_TO_REFIT=1`. Default test runs and tag-CI replays
  do not launch a new full Stan recalibration.
- Deferred to v0.2 and later: correlated hierarchical effects (CHE),
  meta-regression, multivariate outcomes, alternate priors on $`g`$, and
  the empirical-Bayes fallback path. None of these are exposed in v0.1.

### Installation and dependencies

- Imports: `stats`, `splines`, `posterior`, `checkmate`, `rlang`,
  `digest`.
- Suggests: `cmdstanr`, `cli`, `ggplot2`, `ps`, `jsonlite`, `knitr`,
  `rmarkdown`, `testthat`, `metafor`, `withr`, `desc`. CmdStan is
  required to fit models but not to install the package or read cached
  fixtures.
- `Additional_repositories: https://stan-dev.r-universe.dev` so
  `cmdstanr` resolves during installation from a clean R environment.

### Documentation

- The pkgdown site bundles seven applied-track vignettes (A1–A7) and six
  methodological-track vignettes (M1–M6); a five-site cached fit at
  `system.file("examples", "cached_fit_re_smoke.rds", package = "bayesEfron")`
  carries no CmdStan dependency at read time and powers the README Quick
  start and the runnable help examples.
