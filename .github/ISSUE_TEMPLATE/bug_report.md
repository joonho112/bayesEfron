---
name: Bug report
about: Report unexpected behavior in bayesEfron
title: "[bug] "
labels: bug
assignees: ''
---

## Description

A clear, one- or two-paragraph description of what went wrong.

## Reproducible example

A minimal `reprex::reprex()` (or plain R script) that reproduces
the issue. Include:

- The input data shape (`length(theta_hat)`, range of `sigma`, etc.)
- The exact `bayes_efron_fit()` call (or other entry point)
- Any environment variables set
  (`BAYESEFRON_RUN_LIVE`, `BAYESEFRON_NO_GGPLOT2`, etc.)

```r
# paste your reproducible example here
```

## Expected behavior

What you expected to happen.

## Actual behavior

What actually happened. Include the full error message, stack
trace, or unexpected output.

## sessionInfo()

```r
# paste the output of sessionInfo() here
sessionInfo()
```

Include also the output of `cmdstanr::cmdstan_version()` if the
issue involves a Stan fit.
