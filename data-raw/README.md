# Reproducing the example data

Run the scripts from the package root. Neither package installation nor
vignette rendering downloads data or fits a model.

Several scripts take an optional directory, written `/path/to/audit`
below. They write there, as CSV files, tables for checking what they
built.

`krw-discrimination.R` starts from `krw_data.dta` in the discrimination
folder of the public replication package accompanying Walters (2024),
[Empirical Bayes methods in labor
economics](https://sites.google.com/view/christopher-walters/research).
The original experiment is Kline, Rose and Walters (2022), *Systemic
Discrimination Among Large U.S. Employers*, QJE 137(4), 1963-2036; the
97-firm sample follows their (2024) *A Discrimination Report Card*, AER
114(8), 2472-2525,
[doi:10.1257/aer.20230700](https://doi.org/10.1257/aer.20230700).
Install haven to read Stata data. The script uses base R to reproduce OLS
and vacancy-clustered CR1 standard errors and can compare them with the
firm estimates of Lee and Sui (2025),
`data/firm_discrimination_estimates.rds` in their [replication
repository](https://github.com/joonho112/fully-bayes-efron-prior). The
package data are these estimates before any shrinkage.

```sh
Rscript data-raw/krw-discrimination.R /path/to/krw_data.dta \
  /path/to/firm_discrimination_estimates.rds /path/to/audit
```

`twin-towers.R` takes the panel with `I = 0.7` and `R = 9` from
`data/simulated_data_varied_by_I and R.rds` in the [Lee and Sui
replication repository](https://github.com/joonho112/fully-bayes-efron-prior),
commit `f0a57a0235397a420fb5d1562934abbc3f725886`. Its 1,500 true effects
are `deconvolveR::disjointTheta`. Installing deconvolveR enables a check
against those effects. That file contains 15 panels, generated one after
another from the seed 12345. `R` is the ratio of the largest to the
smallest standard error, so the ratio of the variances is its square. The
script copies the panel from the file and does not generate it again,
because newer software need not give the same numbers.

```sh
Rscript data-raw/twin-towers.R "/path/to/simulated_data_varied_by_I and R.rds" \
  /path/to/audit
```

`star.R` builds `star_schools`, the estimates for the 79 kindergarten
schools of Tennessee's Project STAR, from the student file
`STAR_Students.tab` of
[Achilles et al. (2008)](https://doi.org/10.7910/DVN/SIWH9F) on Harvard
Dataverse, which is released under CC0 1.0. Install lme4 and digest. The
script fits one mixed model with a random intercept for each class, checks
the numbers of students, classes and schools, and fails with an error if
the data lines of the file are not the lines that the estimates were built
from. Dataverse generates the header line of the tabular file, and four of
its column names are placeholders, so the script hashes the data lines
apart from the header. If the directory `vignettes/articles/` exists, the
script also writes `class-size-first-stage.rds` there, the same contrasts
under other assumptions about the variance of the classes. If
`data/star_schools.rda` exists, the script compares the new estimates with
it: it fails with an error if they differ by more than 1e-6, and if they
differ by less it writes neither file and prints a message.

```sh
Rscript data-raw/star.R /path/to/STAR_Students.tab /path/to/audit
```

`krw_discrimination`, `twin_towers` and `star_schools` record public
source information and hashes in their `source` attribute. The method and
real-data reanalysis are from
[Lee and Sui (2025)](https://doi.org/10.3390/math13162639).

`raudenbush.R` builds the teacher expectancy data. `grid-comparison.R`
made a comparison of grids for those data with version 0.3.0, which the
tests now read (`tests/testthat/_fixtures/grid_comparison_v030.rds`).
`benchmark/` contains the simulations with known effects; they are not
rerun to render the vignettes.

## The stored fits and the numbers that the documentation reads

`example-fits.R` and `example-evidence.R` need CmdStan, and the second
also needs metafor. Run them in this order from the package root:

```sh
Rscript data-raw/example-fits.R /path/to/audit
Rscript data-raw/example-evidence.R /path/to/audit
```

`example-fits.R` writes the three fits in `data/` with the defaults of
this version: `krw_fit`, `raudenbush_fit` and `twin_towers_fit` (which
sets `M = 8`). Each call runs 4 chains of 1,000 warmup and 3,000 sampling
iterations and stores 1 draw in 15 (`thin = 15`). The script makes every
fit a second time with all its draws, compares the site intervals of the
two, and fails with an error if a fit has not converged, raises a
warning, or has an interval limit more than 2 grid spacings from the same
limit in the full run. It also writes
`inst/examples/twin_towers_comparison.rds`, which compares the fit of the
simulated sites with fits on other ranges, among them the range around
the true effects.

`range-comparison.R` writes `inst/examples/raudenbush_ranges.rds`, the
comparison of three ranges for the teacher expectancy studies that the
examples of two help pages read. It takes a few seconds.

`example-evidence.R` writes `inst/examples/lee_sui_examples.rds`: the
default fit of the callback gaps beside fits with a narrower range and
with 8 spline functions, the rows of the comparison for the simulated
sites, the normal random-effects model fitted by REML to the same
estimates, and the probability of a positive effect for every site. See
`inst/examples/README.md` for what the files contain.

`ebnm-comparison.R` writes `inst/examples/ebnm_comparison.rds`: the stored
fits beside empirical Bayes fits of the same data by the ebnm package, and
beside Efron's estimator for data with equal standard errors, with the
simulated panels of `data-raw/benchmark/`. It needs ebnm and deconvolveR,
makes two fits and takes about two minutes; what it compares was fixed
beforehand in `data-raw/ebnm-protocol.md`. The two fits store 1 draw in
15, as the stored fits do, and so do the fits of
`benchmark/equal-se-run.R`; the fits of `benchmark/panels-run.R` store all
their draws. The script also writes the two tables of the panels in
`inst/verification/`.

## The article about the schools

The article "Small classes in 79 schools" is built for the website only,
from `vignettes/articles/`, and the four files that it reads are in the
same directory and are not part of the installed package. `star.R` writes
`class-size-first-stage.rds`. The two scripts that write the other three
need CmdStan; run them in this order from the package root:

```sh
Rscript data-raw/class-size.R /path/to/audit
Rscript data-raw/benchmark/schools.R 200 14 data-raw/benchmark/results-v040/schools-rows.csv
```

`class-size.R` writes `class-size-fit.rds`, the fit of `star_schools`
with a higher `adapt_delta`, a higher `max_treedepth` and a run twice as
long as the default, and `class-size-evidence.rds`: the same data under
6 seeds with those settings and with the defaults, the fit with the
standard errors that ignore the classes, the fits under the other
assumptions about the variance of the classes, and the normal model
fitted by REML. The script fails with an error unless the stored fit
passes every sampler check except the check of E-BFMI, which is below 0.3
in all 12 runs of these data. `benchmark/schools.R` writes
`class-size-samples.csv`, the table of simulated data sets with the 79
standard errors and known effects; it fits 800 data sets twice each. With
`--summary-only` it makes no fit and computes the table again from the
stored rows.

`README.md` is knitted from `README.Rmd` with `devtools::build_readme()`,
which also makes the figure of the callback gaps. The figure of the
simulated sites is rebuilt with `Rscript data-raw/readme-figures.R` from
the stored fit and the true effects in `twin_towers`.
