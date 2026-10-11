# Compile the Stan programs before the first fit

The first call to
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
compiles the Stan program, which can take up to a minute.
`bayes_efron_compile()` compiles it in advance, so that the first fit
starts sampling immediately. `bayes_efron_clear_cache()` deletes the
compiled models.

## Usage

``` r
bayes_efron_compile(model_family = "RE", quiet = TRUE, force_recompile = FALSE)

bayes_efron_clear_cache(scope = "all")
```

## Arguments

- model_family:

  Model family. Currently only `"RE"` is available.

- quiet:

  If `TRUE` (the default), the compiler output is not shown.

- force_recompile:

  If `TRUE`, compile the programs even if compiled models already exist.

- scope:

  What to delete. `"all"` (the default) deletes the compiled models from
  the cache. The other values are deprecated and are accepted with a
  warning: `"compiled_models"` is equivalent to `"all"`, `"session"`
  deletes no files and only makes the session load the compiled models
  again for the next fit, and `"lock_only"` does nothing.

## Value

`bayes_efron_compile()` invisibly returns the
[`cmdstanr::CmdStanModel`](https://mc-stan.org/cmdstanr/reference/CmdStanModel.html)
object for the current Stan program. `bayes_efron_clear_cache()`
invisibly returns the number of files that it deleted.

## Details

The package has two Stan programs: the current one, and the program of
version 0.3.0, which
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
runs only with `settings = "0.3.0"`. `bayes_efron_compile()` compiles
both.

Compiled models are stored in `tools::R_user_dir("bayesEfron", "cache")`
and reused in later sessions. A model is compiled again after an update
of the package's Stan file or of CmdStan. To store the models somewhere
else, e.g., on a cluster where the home directory is not writable, set
the environment variable `BAYESEFRON_CACHE_ROOT`.

`bayes_efron_compile()` requires the cmdstanr package and a working
installation of CmdStan.

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bayes_efron_compile()

# Delete the compiled models and compile again.
bayes_efron_clear_cache()
bayes_efron_compile()
} # }
```
