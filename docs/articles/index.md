# Articles

### Applied Track

Workflow-oriented walkthroughs from a five-minute landing example
through anatomy of a fit, the
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
bridge, grid choice, diagnostics, plotting, and a worked case study.

- [A1 · Getting
  started](https://joonho112.github.io/bayesEfron/articles/a1-getting-started.md):

  A five-minute first fit on a five-site toy dataset.

- [A2 · Anatomy of a bayesEfron
  fit](https://joonho112.github.io/bayesEfron/articles/a2-anatomy-of-a-fit.md):

  Every field of the bef_fit_re object and every S3 method.

- [A3 · From metafor::escalc() to
  bayesEfron](https://joonho112.github.io/bayesEfron/articles/a3-metafor-bridge.md):

  End-to-end path from group-arm summaries through escalc to a fitted
  bayesEfron model.

- [A4 · Choosing a
  grid](https://joonho112.github.io/bayesEfron/articles/a4-choosing-a-grid.md):

  Side-by-side comparison of the four grid recipes.

- [A5 · Diagnostics in
  practice](https://joonho112.github.io/bayesEfron/articles/a5-diagnostics-in-practice.md):

  How to read diagnose() and decide whether a fit is releasable.

- [A6 · Plotting and
  visualization](https://joonho112.github.io/bayesEfron/articles/a6-plotting-and-visualization.md):

  Four diagnostic-and-summary views of a fitted bayesEfron model.

- [A7 · Case study — Lee–Sui replication
  light](https://joonho112.github.io/bayesEfron/articles/a7-case-study-lee-sui.md):

  One K=50 replication of the Lee–Sui benchmark, end to end.

### Methodological Track

Mathematical background — the empirical-Bayes deconvolution problem, the
Efron log-spline prior, the random-effects hierarchy, grid construction
and spline basis, the Stan model and cache, and the verification ladder
that gates the v0.1 release.

- [M1 · The empirical-Bayes deconvolution
  problem](https://joonho112.github.io/bayesEfron/articles/m1-empirical-bayes-deconvolution.md):

  The Robbins-Efron formulation, prior-art landscape, and motivation for
  the fully Bayesian log-spline approach.

- [M2 · The Efron log-spline
  prior](https://joonho112.github.io/bayesEfron/articles/m2-efron-log-spline-prior.md):

  Derivation, simplex normalisation, half-Cauchy hyperprior, and a
  smoothness sweep.

- [M3 · Bayesian random-effects
  hierarchy](https://joonho112.github.io/bayesEfron/articles/m3-bayesian-hierarchy.md):

  The four-level hierarchy, posterior factorisation, and the package’s
  metadata payload.

- [M4 · Grid construction and spline
  basis](https://joonho112.github.io/bayesEfron/articles/m4-grid-and-basis.md):

  Mathematical specification of the four grid recipes and the
  natural-cubic-spline basis.

- [M5 · Stan model, sampler, and
  cache](https://joonho112.github.io/bayesEfron/articles/m5-stan-and-cache.md):

  The locked Stan model, the seven-field data block, and the two-tier
  compile cache.

- [M6 · Verification and
  calibration](https://joonho112.github.io/bayesEfron/articles/m6-verification-and-calibration.md):

  The verification tier ladder and the Lee-Sui aggregate coverage
  evidence that gates v0.1.0.
