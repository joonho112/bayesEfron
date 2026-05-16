# Pre-compile the bayesEfron Stan model

Pre-warm the CmdStan compilation cache used by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).
The function compiles (or reattaches a cached build of) the locked v0.1
random-effects Stan model, stores the resulting
[`cmdstanr::CmdStanModel`](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)
in the in-session cache, and runs a tiny post-compile smoke check to
confirm the binary is callable.

Calling `bayes_efron_compile()` before the first
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
shifts the (potentially long) compilation cost out of the fit pipeline,
which is useful when fitting interactively or under a wall-clock budget.
It is otherwise optional:
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
triggers the same cache mechanism on first use.

## Usage

``` r
bayes_efron_compile(
  model_family = "RE",
  quiet = TRUE,
  force_recompile = FALSE,
  seed_for_check = 42L
)
```

## Arguments

- model_family:

  Character scalar. v0.1 supports `"RE"` only.

- quiet:

  Logical. If `TRUE` (default), compile and smoke-check output is
  suppressed during cache warming.

- force_recompile:

  Logical. If `TRUE`, bypasses cached artifacts and recompiles from Stan
  source. Defaults to `FALSE`.

- seed_for_check:

  Non-negative integer scalar used for the synthetic post-compile smoke
  check. Defaults to `42L`.

## Value

Invisibly, a
[`cmdstanr::CmdStanModel`](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)
reference attached to the cache entry. The return value is rarely used
directly; it is returned to support advanced workflows that want to
drive the model object outside the package's pipeline.

## Details

The cache lives at the location given by the environment variable
`BAYESEFRON_CACHE_ROOT` (with a sensible per-user default if not set).
The lookup key combines the Stan source SHA-256, the package version,
the CmdStan version, and the platform, so a cached binary is only reused
when all of those match. The post-compile smoke check uses the same
fixed internal sampler initialization as
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
(`init = 0.5`).

This function requires the `cmdstanr` package and a working CmdStan
toolchain. Both are listed in `Suggests:` rather than `Imports:` so the
package can be installed and documented without them.

## See also

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  for the user-facing fit pipeline.

- [`bayes_efron_clear_cache()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_clear_cache.md)
  for cache maintenance and stale lock recovery.

## Examples

``` r
if (FALSE) { # \dontrun{
# Pre-warm the cache so the next fit skips the compile cost.
bayes_efron_compile()
} # }
```
