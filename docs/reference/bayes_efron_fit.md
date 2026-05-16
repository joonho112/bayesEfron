# Fit the Bayesian Efron random-effects model

Fit the fully Bayesian Efron log-spline prior to a univariate
random-effects meta-analytic deconvolution problem. Given per-site
effect estimates and their within-study standard errors,
`bayes_efron_fit()` returns posterior site-effect summaries together
with a continuous estimate of the underlying mixing distribution.

The function targets applied meta-analysts who already have point
estimates and standard errors on a comparable scale (for example, the
`yi` and `sqrt(vi)` columns of an
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
object). One call performs the full eight-stage pipeline: input
validation, grid construction, Stan-data preparation, cache-backed model
retrieval, NUTS sampling, draw extraction, generated-quantity
postprocessing, and assembly into a validated `bef_fit_re` object.

Computation is delegated to CmdStan through `cmdstanr`; the first call
in a session typically spends most of its time compiling the Stan
program. Subsequent calls reuse the on-disk and in-session caches as
long as the model source, package version, CmdStan version, and platform
match.

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
  iter_warmup = 1000L,
  iter_sampling = 3000L,
  adapt_delta = 0.9,
  seed = NULL,
  keep_cmdstan_fit = FALSE
)
```

## Arguments

- theta_hat:

  Numeric vector of per-site effect estimates on a common scale (e.g.
  mean differences, log odds ratios).

- sigma:

  Numeric vector of strictly positive per-site standard errors, on the
  same scale as `theta_hat` and the same length.

- ...:

  Reserved for future expansion; must be empty in v0.1.

- grid_method:

  Character grid recipe. One of `"paper_realdata"` (default),
  `"paper_simulation"`, `"paper_sensitivity"`, or
  `"kl_target_experimental"`. See **Details**.

- L:

  Integer grid length (number of discrete support points). Defaults to
  `101L`. The package's verification ledger is calibrated at this
  default.

- expansion:

  Numeric, non-negative. Range-relative expansion factor that widens the
  grid endpoints beyond the observed range of `theta_hat`. Defaults to
  `0.5` (50 percent expansion).

- M:

  Integer natural-cubic-spline degrees of freedom. Defaults to `6L`. The
  verification ledger is calibrated at this default.

- theta_true:

  Numeric oracle vector of latent site effects, required by
  `"paper_simulation"` and `"paper_sensitivity"` and ignored by the
  other recipes. Same length as `theta_hat`.

- bound_expansion:

  Numeric, oracle-bound expansion factor used only by
  `"paper_sensitivity"`. `NULL` falls back to the recipe default of
  `0.5`.

- model_family:

  Character scalar. v0.1 supports `"RE"` only. Other model families are
  deferred to v0.2+ per the package blueprint.

- chains:

  Integer number of MCMC chains. Defaults to `4L`.

- iter_warmup:

  Integer warmup iterations per chain. Defaults to `1000L`.

- iter_sampling:

  Integer post-warmup iterations per chain. Defaults to `3000L`.

- adapt_delta:

  NUTS target acceptance statistic in `(0, 1)`. Defaults to `0.9`.

- seed:

  Integer seed or `NULL`. If `NULL`, an integer seed is auto-generated
  and recorded on `fit$metadata$seed`.

- keep_cmdstan_fit:

  Logical. If `TRUE`, the raw
  [`cmdstanr::CmdStanMCMC`](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)
  handle is retained at `fit$cmdstan_fit` for advanced use. Defaults to
  `FALSE` so the returned object is small enough to save and share.

## Value

An S3 object of class `c("bef_fit_re", "bef_fit")` with three top-level
fields:

- `draws` — a
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
  of MCMC draws for the model parameters and generated quantities.

- `metadata` — a named list with exactly 13 fields:

  - `model_family` — `"RE"` for v0.1.

  - `grid_method` — the recipe used.

  - `seed` — effective integer seed (auto-generated if not supplied).

  - `cmdstan_version` — CmdStan version string.

  - `stan_file_sha256` — SHA-256 of the locked Stan source.

  - `data_list` — the seven-field Stan data block sent to CmdStan.

  - `runtime_seconds` — sampler wall-clock.

  - `mean_g_summary`, `var_g_summary` — posterior summaries of
    functionals of the mixing distribution \\g\\.

  - `theta_summary` — posterior summaries of latent site effects
    \\\theta_i\\.

  - `theta_rep_draws` — replicated effects for posterior predictive
    checks.

  - `effective_params_summary`, `log_marginal_likelihood_summary` —
    model-quality summaries.

- `posterior` — a tidy posterior representation used by the S3 methods.

Four additional durable payloads are stored as attributes of
`fit$metadata`: `sd_g_summary`, `diagnostics`, `diagnostic_skipped`, and
`sampler_diagnostics_failed`. Access them through
[`summary.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
and
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
rather than directly.

When `keep_cmdstan_fit = TRUE`, `fit$cmdstan_fit` carries the raw
[`cmdstanr::CmdStanMCMC`](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)
handle.

## The fitted four-level hierarchy

For sites \\i = 1, \ldots, K\\, \$\$\hat\theta_i \mid \theta_i \sim
\mathcal{N}(\theta_i, \sigma_i^2),\$\$ \$\$\theta_i \mid g
\overset{\text{iid}}{\sim} g,\$\$ where \\g\\ is the mixing distribution
of latent site effects. The package places a log-spline prior on the
discretized \\g\\ with a half-Cauchy hyperprior on the smoothness
precision \\\lambda\\; full derivations are in the methodological
vignettes.

## Grid recipes

The four `grid_method` choices control how the discrete support of \\g\\
is constructed:

|  |  |  |
|----|----|----|
| Recipe | Needs `theta_true`? | Use when |
| `"paper_realdata"` | No | Real-data analysis with no oracle. |
| `"paper_simulation"` | Yes | Simulation with known truth, matched paper. |
| `"paper_sensitivity"` | Yes | Bound-expansion sensitivity, paper rule. |
| `"kl_target_experimental"` | No | KL-target tuning (experimental, heteroscedastic). |

The two oracle-requiring recipes refuse to run without a numeric
`theta_true`. The experimental recipe emits a once-per-session
disclaimer about its KL calibration.

## Sampler defaults

v0.1 fixes the CmdStan initialization at `init = 0.5` to keep release
fits reproducible; this is **not** a user-facing tuning argument.
`parallel_chains` defaults to `chains` but can be overridden by setting
the environment variable `BAYESEFRON_PARALLEL_CHAINS` (or the unprefixed
`PARALLEL_CHAINS`) to a positive integer that does not exceed `chains`.

## Reproducibility

If `seed` is `NULL`, the function auto-generates an integer seed from
the current time and stores it on `fit$metadata$seed` so the fit can be
re-played later. Pass an explicit `seed` for fully reproducible runs.

## The 13-field metadata contract

`fit$metadata` always has exactly the 13 named fields listed in
**Value** below. Four additional payloads are stored as attributes of
`fit$metadata` and surfaced through
[`summary.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md)
and
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md):
`sd_g_summary`, `diagnostics`, `diagnostic_skipped`, and
`sampler_diagnostics_failed`.

## See also

- [`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
  for converting
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  objects and plain lists into the package input class.

- [`make_efron_grid()`](https://joonho112.github.io/bayesEfron/reference/make_efron_grid.md)
  for building grids outside the fitting pipeline (for example, to
  inspect a recipe before fitting).

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  for pre-warming the Stan model cache.

- [`summary.bef_fit()`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md),
  [`confint.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/bayesEfron-methods.md),
  [`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md),
  [`plot.bef_fit_re()`](https://joonho112.github.io/bayesEfron/reference/plot.bef_fit_re.md)
  for inspecting the returned object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Five-site smoke fit. Requires a working CmdStan installation.
theta_hat <- c(-0.21, 0.04, 0.19, 0.38, 0.61)
sigma     <- c( 0.18, 0.15, 0.22, 0.19, 0.24)

fit <- bayes_efron_fit(
  theta_hat     = theta_hat,
  sigma         = sigma,
  L             = 51L,
  M             = 3L,
  chains        = 1L,
  iter_warmup   = 150L,
  iter_sampling = 4L,
  seed          = 1234L
)

summary(fit)
confint(fit, type = "theta")
diagnose(fit)
plot(fit, type = "caterpillar")
} # }
```
