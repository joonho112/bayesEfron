# Compile the Stan model ahead of time, or remove the compiled model

The first call to
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
compiles the Stan model, which can take up to a minute.
`bayes_efron_compile()` does that step on its own, so that a later fit
starts sampling at once. `bayes_efron_clear_cache()` deletes the
compiled model.

## Usage

``` r
bayes_efron_compile(model_family = "RE", quiet = TRUE, force_recompile = FALSE)

bayes_efron_clear_cache(scope = "all")
```

## Arguments

- model_family:

  Character string. `"RE"` is the only model available.

- quiet:

  Logical. If `TRUE` (the default), the compiler output is not shown.

- force_recompile:

  Logical. If `TRUE`, the model is compiled even when a compiled copy
  exists. Defaults to `FALSE`.

- scope:

  Character string. `"all"` (the default) deletes the cached files. The
  values `"session"`, `"compiled_models"` and `"lock_only"` are accepted
  for scripts written for earlier versions and give a warning;
  `"lock_only"` no longer does anything.

## Value

`bayes_efron_compile()` returns the
[`cmdstanr::CmdStanModel`](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)
object, invisibly. `bayes_efron_clear_cache()` returns the number of
files it deleted, invisibly.

## Details

The compiled model is stored in
`tools::R_user_dir("bayesEfron", "cache")` and reused in later sessions.
The file name includes a hash of the Stan file, the CmdStan version and
the machine architecture, so the model is compiled again after an update
of the package's Stan file or of CmdStan. Set the environment variable
`BAYESEFRON_CACHE_ROOT` to store it somewhere else, for example on a
cluster where the home directory is not writable.

`bayes_efron_compile()` needs the cmdstanr package and a working CmdStan
installation.

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bayes_efron_compile()

# Start again from the Stan file.
bayes_efron_clear_cache()
bayes_efron_compile()
} # }
```
