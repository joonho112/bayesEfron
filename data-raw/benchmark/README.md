# Benchmark with known effects

The scripts and data behind the two tables in `inst/verification/`, which
`vignette("implementation")` shows. Run the scripts from the package root;
they need CmdStan and the packages devtools and digest.

- `benchmark-coverage.R` fits the model to 20 data sets for each of five
  numbers of sites (50, 100, 200, 500 and 1,500) and computes the coverage
  of the 90% intervals with the convergence diagnostics. It writes
  `benchmark-coverage.csv`, the table that the package ships, and
  `benchmark-per-replication.csv`, with one row for each fit. The 100 fits
  take one to two hours.
- `benchmark-remedies.R` refits the data sets with 1,500 sites whose fits
  did not converge, with a narrower grid and with a larger maximum tree
  depth, and writes `benchmark-remedies.csv`. It reads the per-replication
  results of the first script.

Both write to `data-raw/benchmark/results/` unless another directory is
given as the first argument.

The data are in `lee_sui_K100.rds`, `lee_sui_K200.rds`, `lee_sui_K500.rds`
and `lee_sui_K1500.rds`; the file for 50 sites is
`tests/testthat/_fixtures/lee_sui_K50.rds`, because a test reads it. Each
file holds 20 data sets with the estimates, their standard errors, the true
effects and the grid settings of Lee and Sui (2025, Appendix E). The data
sets with 50 to 500 sites are those of the analysis in that appendix:
stratified samples of the sites of one simulated panel of 1,500. The 20
data sets with 1,500 sites were simulated with the data-generating rule of
the paper's replication code (`I = 0.7` and `R = 9`, so that the largest
standard error is 9 times the smallest): they share the true effects, and
each has its own assignment of the standard errors to the sites and its own
errors.

The element `metadata` of each file names the files from which the data
were assembled and gives their hashes. Those files are not in this
repository, and the script that assembled the data is not distributed; the
two scripts here start from the data files. The `seed` in a data file is
the seed with which that data set was generated. The sampler has its own:
both scripts fit replication `r` with the seed `20260509 + r`.
