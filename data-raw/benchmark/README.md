# Simulated data sets with known effects

The scripts and data behind the tables in `inst/verification/`, which
`vignette("implementation")` sets side by side, and behind the table of
simulated schools in the article "Small classes in 79 schools". Run the
scripts from the package root; they need CmdStan and the package devtools.

## The data of Lee and Sui (2025, Appendix E)

`lee_sui_K100.rds`, `lee_sui_K200.rds`, `lee_sui_K500.rds` and
`lee_sui_K1500.rds` each contain 20 data sets with the estimates, their
standard errors, the true effects and the grid settings of that appendix;
the file for 50 sites is `tests/testthat/_fixtures/lee_sui_K50.rds`,
because a test reads it. The data sets with 50 to 500 sites are the data
sets of the analysis in the appendix: stratified samples of the sites of
one simulated panel of 1,500. The 20 data sets with 1,500 sites were
simulated with the data-generating process of the paper's replication code
(`I = 0.7` and `R = 9`, so that the largest standard error is 9 times the
smallest): they have the same true effects, and each has its own
assignment of the standard errors to the sites and its own errors.

The element `metadata` of each file names the files from which the data
were assembled and gives their hashes. Those files are not in this
repository, and the script that assembled the data is not distributed.
The `seed` in a data file is the seed with which that data set was
generated.

Version 0.3.0 shipped two tables computed from these data, the coverage of
the intervals by number of sites and the refits of the large data sets
that did not converge. They described the settings of that version and
are not part of the package any more; the scripts that made them are in
the source of version 0.3.0.

## Independent panels

The 20 data sets with 1,500 sites above have one set of true effects in
common, so they are not independent replications. `panels-make.R`
generates panels that each have their own true effects, drawn from the two
uniform populations of the Twin Towers example, with the same reliability
and the same ratio of the standard errors:

```sh
Rscript data-raw/benchmark/panels-make.R            # panels_K1500.rds, 20 panels
Rscript data-raw/benchmark/panels-make.R 100 100    # panels_K100.rds, 100 panels
```

`panels-run.R <arm> [outdir] [K]` fits them with the defaults of this
version (`default`), with 8 spline functions (`M8`) or with the settings
of version 0.3.0 (`v030`), and stores for every fit a row of diagnostics,
cost and recovery, and the estimate of g with the table of the sites.
`equal-se-run.R <M>` adds errors with a standard deviation of 1 to the
true effects of the large panels and compares the fit with Efron's
estimator. A run that was interrupted continues where it left off.

`results-v040/` contains what these runs wrote for version 0.4.0 on
2026-10-09: the 20 large panels in the three arms, the 100 small panels
with the defaults, and the runs with equal standard errors. From them
`data-raw/ebnm-comparison.R` builds the tables
`inst/verification/panels-convergence.csv` and
`inst/verification/panels-methods.csv`. The fits of the three arms of the
large panels took 26, 19 and 53 minutes in all (the sums of
`elapsed_seconds`); the times in the files were measured on a machine that
was running other jobs at the same time.

## Nineteen studies

`few-studies.R` simulates data sets with the standard errors of the 19
teacher expectancy studies and true effects from three distributions, fits
each with the defaults and with `adapt_delta = 0.99`, and writes
`inst/verification/few-studies.csv`, the table of
`vignette("teacher-expectancy")`. It makes 1,200 fits and takes about ten
minutes on 14 cores.

`results-v040/few-studies-rows.csv` contains the results of all 600 data
sets, four analyses for each. The rows called `known g` are computed
without a fit. For the third distribution an effect is a grid point moved
uniformly within its cell, so those rows use the prior distribution that
is uniform within each cell: `histogram-posterior.R` has the formulas, and
`histogram-posterior-check.R` compares them with numerical integration and
computes the 200 rows again from their seeds. This command computes those
200 rows and their row of the table again without a fit, and leaves every
other row, `generated_on` and `cmdstan_version` as they are:

```sh
Rscript data-raw/benchmark/few-studies.R --analytic-only
```

A second and a third argument give the rows file and the table.

## Seventy-nine schools

`schools.R` simulates data sets with the standard errors of the 79 schools
of `star_schools` and known effects, fits each of them twice, and writes
`vignettes/articles/class-size-samples.csv`, the table of the article
"Small classes in 79 schools". `results-v040/schools-rows.csv` contains its
rows for every data set. `data-raw/README.md` gives the command.
