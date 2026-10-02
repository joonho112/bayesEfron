# Benchmark coverage

`benchmark-coverage.csv` holds the coverage of the 90% intervals on the
benchmark of Lee and Sui (2025), in which the true effects are known. It is
the table shown in `vignette("implementation")`, which also describes the
design.

Each row is a number of sites `K`. `coverage` is the share of sites whose
true effect lies inside the interval, over `n_replications` data sets, and
`mean_width` is the mean width of the intervals. `max_rhat` and
`min_ess_bulk` are the worst values over the fits, `n_rhat_above_1.05` is
the number of fits whose largest R-hat exceeded 1.05, and
`total_divergences` and `total_treedepth_hits` are totals over the fits.
`generated_on` and `cmdstan_version` say when and with which CmdStan the
table was computed.

`coverage_se` is the standard deviation of the coverage between the data
sets of one size, divided by the square root of their number. It measures
the variation between those data sets and nothing more. The data sets of
50 to 500 sites are samples of the sites of one simulated panel, and those
of 1,500 sites share one set of true effects and differ in the errors and
in the assignment of the standard errors to the sites. The variation that
another panel or another set of true effects would bring is not in it.

`benchmark-remedies.csv` holds, for the fits with 1,500 sites whose largest
R-hat exceeded 1.05, the same diagnostics and the coverage with the default
settings, with a narrower grid (`expansion = 0.05`) and with a maximum tree
depth of 12. `replication` is the number of the data set.

The tables are not recomputed when the package is checked; computing them
takes more than an hour. The scripts that compute them and the data for 100
to 1,500 sites are in `data-raw/benchmark/` in the source repository. The
data for 50 sites are in `tests/testthat/_fixtures/`.

The element `metadata` of each data file names the files from which the
data were assembled and gives their hashes. Those files are not part of
the package, and the script that assembled the data is not distributed;
the benchmark scripts start from the data files. The `seed` in a data file
is the seed with which that data set was generated. The sampler has its
own: the fit of replication `r` uses `20260509 + r`.
