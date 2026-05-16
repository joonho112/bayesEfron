# Clear bayesEfron compilation cache artifacts

Remove selected in-session and on-disk cache artifacts used by the
CmdStan compilation cache. Intended for cache maintenance, stale-lock
recovery after an interrupted compile, and forced rebuilds during
development.

Choose the smallest scope that resolves the issue at hand. The
`"lock_only"` scope is the safest and is appropriate for clearing a
stale lock left behind by a killed process. `"all"` is the largest scope
and will force every subsequent fit in this session (and on this
machine, until the cache is repopulated) to recompile from scratch.

## Usage

``` r
bayes_efron_clear_cache(
  scope = c("lock_only", "session", "compiled_models", "all")
)
```

## Arguments

- scope:

  Character scalar. One of `"lock_only"` (default), `"session"`,
  `"compiled_models"`, or `"all"`.

## Value

Invisibly, a named integer vector with elements `lock_files`,
`session_keys`, `disk_models`, and `disk_sidecars`. Each element counts
the number of artifacts removed; the contract is fixed even when the
scope leaves a given counter at zero.

## Scopes

|  |  |  |
|----|----|----|
| Scope | Removes | Use when |
| `"lock_only"` | Stale lock file | Resuming after an interrupted compile that left a lock behind. |
| `"session"` | The in-session cache only | Forcing the current session to re-attach to disk-cached binaries. |
| `"compiled_models"` | Cache-entry binaries, sidecar JSON, and companion Stan-source copies under the current cache format directory; preserves the cache directory and lock file | Forcing a recompile while keeping the cache root intact. |
| `"all"` | The entire cache root, plus the in-session cache | Resetting the cache to a clean slate. |

This function does **not** acquire the cache lock. Avoid calling the
disk-clearing scopes while another R process is compiling a bayesEfron
model in the same cache root.

The cache root location is controlled by the environment variable
`BAYESEFRON_CACHE_ROOT` (with a sensible per-user default if unset).

## See also

- [`bayes_efron_compile()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_compile.md)
  for repopulating the cache after a clear.

- [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
  for the user-facing fit pipeline that consumes the cache.

## Examples

``` r
# Non-destructive: clear a stale lock if one is present.
bayes_efron_clear_cache("lock_only")
#> bayesEfron cache clear (lock_only): lock_files=0, session_keys=0, disk_models=0, disk_sidecars=0.

if (FALSE) { # \dontrun{
# Destructive: drop the current session's cache entries.
bayes_efron_clear_cache("session")

# Most destructive: reset the entire cache root.
bayes_efron_clear_cache("all")
} # }
```
