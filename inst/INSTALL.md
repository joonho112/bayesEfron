# Installing bayesEfron

This guide separates R package installation from CmdStan setup. You can install
and load `bayesEfron`, read help files, build grids, and inspect input classes
without CmdStan. Calls to `bayes_efron_fit()` and `bayes_efron_compile()` need
`cmdstanr`, a C++ toolchain, and CmdStan.

## Requirements

- R 4.1.0 or newer.
- Package dependencies listed in `DESCRIPTION`.
- A working C++ toolchain for CmdStan.
- CmdStan 2.34.0 or newer for fitting and live verification.

## Install The R Package

From GitHub:

```r
install.packages("remotes")
remotes::install_github("joonho112/bayesEfron")
```

With `pak`:

```r
install.packages("pak")
pak::pkg_install("joonho112/bayesEfron")
```

From a local source checkout:

```r
install.packages(
  "/absolute/path/to/bayesEfron-R-package",
  repos = NULL,
  type = "source"
)
```

## Install cmdstanr And CmdStan

Install `cmdstanr` from the Stan R-universe repository:

```r
install.packages(
  "cmdstanr",
  repos = c(
    "https://stan-dev.r-universe.dev",
    "https://cloud.r-project.org"
  )
)
```

Check the toolchain and install CmdStan:

```r
cmdstanr::check_cmdstan_toolchain(fix = FALSE)
cmdstanr::install_cmdstan(cores = 2)
cmdstanr::cmdstan_version()
```

If the toolchain check fails, follow the platform-specific instructions from
the CmdStan and cmdstanr documentation, then rerun the checks above.

## Prewarm The Stan Cache

The first fit compiles the locked random-effects Stan model. To compile before
an analysis session:

```r
bayesEfron::bayes_efron_compile(quiet = FALSE)
```

The cache key records the model source digest, package version, CmdStan version,
R platform, and cache format version. A later fit can reattach the cached model
when these fields match.

## Cache Location

By default, `bayesEfron` chooses a user cache directory through R's standard
cache conventions. To make cache behavior explicit for a project or CI job, set:

```sh
export BAYESEFRON_CACHE_ROOT=/absolute/path/to/bayesEfron-cache
```

For a temporary CI cache, point this variable at the CI workspace cache
directory and restore it between jobs.

## Troubleshooting

If a previous compile was interrupted and left a stale lock:

```r
bayesEfron::bayes_efron_clear_cache("lock_only")
```

If a compiled binary is stale or incompatible with the current toolchain:

```r
bayesEfron::bayes_efron_clear_cache("compiled_models")
bayesEfron::bayes_efron_compile(quiet = FALSE, force_recompile = TRUE)
```

To reset all bayesEfron cache artifacts for the active cache root:

```r
bayesEfron::bayes_efron_clear_cache("all")
```

Avoid clearing compiled model files while another R process is compiling the
same model.

## Live Test Gates

Default checks skip live CmdStan work. Opt in when the local machine has
CmdStan configured:

```sh
BAYESEFRON_RUN_LIVE=1 Rscript -e 'testthat::test_local(".")'
BAYESEFRON_RUN_FULL_LIVE=1 Rscript -e 'testthat::test_local(".")'
BAYESEFRON_RUN_PARITY=1 Rscript -e 'testthat::test_local(".")'
```

For CI gating, use `R CMD check` or the existing GitHub Actions
`r-lib/actions/check-r-package@v2` workflow as the pass/fail authority. Do not
add CI gates that rely only on the process exit status from bare
`Rscript -e 'devtools::test(...)'`; prior negative-control work observed that
this can return exit status 0 even when test failures are reported. If a custom
scripted test gate is needed, wrap the test call so failures are converted into
an explicit non-zero exit.

The full live gate is intended for nightly, release-candidate, or manually
supervised verification because it can be slow.

Fresh Tier 3 full-live refits are intentionally double-gated. When
`BAYESEFRON_RUN_FULL_LIVE=1` is set without
`BAYESEFRON_TIER3_FULL_LIVE_MATRIX`, the test harness requires
`BAYESEFRON_TIER3_OK_TO_REFIT=1` before launching new K50-K1500 Stan fits. This
does not affect replay checks that provide an accepted matrix.

At v0.1.0, the Tier 3 Lee-Sui fixtures are present as locked
all-20-replication bundles. Setting `BAYESEFRON_RUN_FULL_LIVE=1` opts into
supervised full-live work over K50, K100, K200, K500, and K1500 with all 20
replications per K. The verification ledger records those five Tier 3
aggregate coverage rows as active.

For regression replay against a previously computed coverage matrix, set
`BAYESEFRON_TIER3_FULL_LIVE_MATRIX` to a CSV path with the schema
documented in `tests/testthat/test-tier3-live-harness.R`. The matrix CSV is
not shipped with the package; supply one produced from a prior run or use
the fresh-refit path (requires `BAYESEFRON_TIER3_OK_TO_REFIT=1` in
addition to `BAYESEFRON_RUN_FULL_LIVE=1`).

The accepted Tier 3 evidence is benchmark-scoped release evidence for
`theta_rep` interval calibration on the archived Lee-Sui all-20-replication
fixtures under the v0.1 default configuration. Extending the same empirical
coverage claim to new grids, spline degrees, data regimes, or shortened
sampling settings requires re-verification.

The repaired K1500 acceptance was close to the lower bound: aggregate
`theta_rep` interval coverage was 0.876433333333333 in the pre-specified
`[0.87, 0.92]` interval over 20 completed replications. Treat this as accepted
release evidence, not as a large-margin robustness claim.

The repository workflow `.github/workflows/check.yml` uses three envelopes:

- default push and pull-request checks leave live gates unset;
- pull requests labelled `live-ci` and manual dispatches with `live = smoke`
  run the Ubuntu release live-smoke envelope;
- version tags, scheduled runs, and manual dispatches with `live = full` run
  the Ubuntu release full-live envelope; the workflow prepares a CI-local
  replay matrix when one is available in the secrets store and otherwise
  triggers the fresh-refit path.

### Tier 3 Reproducibility Caveat

The Tier 3 `lee_sui_K*.rds` fixtures are canonical to the archived Lee-Sui
appendix object `small-k-scenarios_simulation_results.rds`, not to a fresh replay
of the public `Part_06_Small-K Scenarios Analysis.R` script alone. Each fixture
stores all 20 archived replications for its K, plus replication-1 top-level
aliases for smoke-mode compatibility. The public Part_06 script is retained in
fixture metadata as public replication context, but it uses a median split of
`theta_true` and sorted sampled indices. The archived appendix object was
produced by the original development script using a `theta_true < 0` versus
`theta_true >= 0` split and unsorted sampled indices. The package generator
follows the archived appendix object so that the fixture vectors match the
published small-K panels byte-for-byte.

### Tier 3 Sampler Diagnostic Caveat

The accepted full-live Tier 3 evidence completed the release-blocking
all-20-replication matrix and the repaired K1500 full-settings run completed
all 20 replications with zero divergences. The accepted replication fits still
recorded sampler diagnostic flags for `rhat`, `ess_bulk`, and `ess_tail`. These
flags are retained as diagnostic metadata for the release-calibration tier. They
do not invalidate the accepted aggregate coverage rows, but future fresh
full-live runs should continue to surface and review them.
