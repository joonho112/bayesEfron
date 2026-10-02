# Installing bayesEfron and CmdStan

bayesEfron is installed from GitHub.

```r
# install.packages("remotes")
remotes::install_github("joonho112/bayesEfron")
```

`install_github()` does not build the vignettes unless it is asked to.
They can be read on the
[package website](https://joonho112.github.io/bayesEfron/). To have them in
R as well, install knitr and rmarkdown, make sure that Pandoc is available
(RStudio includes it; otherwise see <https://pandoc.org/installing.html>),
and install with

```r
remotes::install_github("joonho112/bayesEfron", build_vignettes = TRUE)
```

The package fits its model with CmdStan, through the cmdstanr package.
Neither is on CRAN. The examples and the vignettes run without them,
because they use stored results; they are needed to fit a model.

```r
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

CmdStan is built from source, so a C++ toolchain is needed first.

- macOS: the Xcode command line tools, installed with `xcode-select --install`.
- Windows: Rtools; `cmdstanr::check_cmdstan_toolchain(fix = TRUE)` sets up
  what is missing.
- Linux: `g++` and `make`, for example from the package `build-essential`.

The cmdstanr documentation,
<https://mc-stan.org/cmdstanr/articles/cmdstanr.html>, has the details.

To check the installation and compile the model of bayesEfron ahead of the
first fit:

```r
cmdstanr::cmdstan_version()
bayesEfron::bayes_efron_compile()
```

The compiled model is kept in `tools::R_user_dir("bayesEfron", "cache")`.
Where the home directory is not writable, as on some clusters, set the
environment variable `BAYESEFRON_CACHE_ROOT` to a directory that is.
