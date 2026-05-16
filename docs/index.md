# bayesEfron

**Fully Bayesian Efron log-spline deconvolution for univariate
random-effects meta-analysis with heteroscedastic standard errors.**

You arrive with per-site effect estimates $`\hat\theta_i`$ and
within-study standard errors $`\sigma_i`$, one pair per study, with
$`\sigma_i`$ free to vary across sites. The estimates are noisy
realizations of latent site effects $`\theta_i`$ drawn from an unknown
mixing distribution $`g(\theta)`$. The conventional random-effects
meta-analysis collapses $`g`$ to a Gaussian summarized by a single
heterogeneity parameter, which discards information whenever the latent
effects are skewed, heavy-tailed, or multimodal.

`bayesEfron` returns two posterior objects from a single fit. The first
is a continuous nonparametric estimate of $`g(\theta)`$ itself,
expressed through Efron’s log-spline prior on a fixed grid of
effect-size values. The second is the per-site posterior
$`\theta_i \mid \hat\theta_i, \sigma_i`$ — empirical-Bayes shrunken
estimates that borrow strength under the estimated $`g`$. Deconvolution
is the operation that recovers $`g`$ from the noisy
$`(\hat\theta_i, \sigma_i)`$ pairs; it is what makes the per-site
posteriors honest about the population they are drawn from.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("joonho112/bayesEfron")
```

CmdStan is required to fit models but not to install the package or to
load cached fixtures. See
`system.file("INSTALL.md", package = "bayesEfron")` for the CmdStan
setup guide.

## Quick start

``` r

library(bayesEfron)
fit <- readRDS(system.file("examples", "cached_fit_re_smoke.rds", package = "bayesEfron"))
summary(fit)
confint(fit, type = "theta")
plot(fit, type = "caterpillar")
```

The cached fixture is a five-site smoke fit produced by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
with `L = 51`, `M = 3`, one chain, and 100 sampling iterations; it
carries no CmdStan dependency at read time. The verbatim `summary(fit)`
output:

    <summary.bef_fit>

    Prior g:
      mean: 0.174
      var:  0.2094
      sd:   0.4549

    Diagnostics:
      Rhat:              1.203
      ESS bulk:          12.53
      ESS tail:          105
      Divergences:       0
      Max treedepth:     0
      Effective params:  mean 0.3754, sd 0.3722
      Log marginal lik.: mean -2.433, sd 0.443
      Runtime:           0.2906 sec
      Stan SHA-256:      57d1f55c8ccd721800193e61eb247f315be3afd401a19f58a51f84c103ceca3e

    Theta summary:
      site      mean        sd     lower     upper       map
      1      -0.1948    0.1713    -0.456     0.036   -0.1956
      2       0.0401    0.1471   -0.1608       0.2   0.03862
      3       0.1808    0.2128   -0.1952    0.4968    0.1787
      4       0.3652    0.1869   0.06716    0.6608     0.363
      5       0.5658    0.2197    0.2984    0.9216    0.5847

The **Prior g** block holds the posterior-mean summary of the mixing
distribution. **Theta summary** is the per-site posterior table.
**Diagnostics** are at smoke-scale by construction; production fits use
the package defaults (`L = 101`, `M = 6`, four chains, 1000 warmup, 3000
sampling).

## Key features

- **Fully Bayesian inference.** The spline coefficients for
  $`g(\theta)`$ are sampled jointly with the latent $`\theta_i`$,
  propagating uncertainty about $`g`$ into every per-site posterior —
  not a plug-in empirical-Bayes point estimate of $`g`$ followed by a
  separate conditional step.
- **Exact heteroscedasticity.** Each site enters the likelihood with its
  own $`\sigma_i`$; the model never collapses the within-study standard
  errors to a single pooled value or a meta-analytic typical variance.
- **Two-tier compile cache.** A first fit compiles the Stan model and
  stores the binary in both an in-session cache and an on-disk cache
  keyed on model source, package version, CmdStan version, and platform.
  Later fits with matching keys skip compilation entirely.

## Where to next

- **A1 · Getting started** — a runnable five-minute walkthrough of the
  same five-site fixture, with both inferential targets read off the fit
  object.
- **A3 · From
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  to bayesEfron** — the end-to-end aggregate-data path, starting from
  group treatment and control summaries through a working fit and its
  outputs.
- **M1 · The empirical-Bayes deconvolution problem** — the formal
  motivation for the log-spline prior, the role of fully Bayesian
  inference over Efron’s original empirical-Bayes formulation, and the
  connection to the broader literature.

## Citation

To cite `bayesEfron`, see `citation("bayesEfron")` for the canonical
reference.

## Funding and disclaimer

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama. The opinions expressed are those of the authors and do not
represent views of the Institute or the U.S. Department of Education.
