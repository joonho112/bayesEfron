# Test fixtures

`grid/` holds the inputs and the expected grids for the three grid methods
that follow Lee and Sui (2025). The tests rebuild each grid with
`make_efron_grid()` and compare it to the stored one.

`deconvolveR_baseline.rds` holds empirical Bayes estimates computed with
`deconvolveR::deconv()` for three simulated data sets (K = 20, 50, 100).

`lee_sui_K50.rds` holds the 20 data sets with 50 sites of the benchmark of
Lee and Sui (2025, Appendix E). The files for 100, 200, 500 and 1,500 sites
and the scripts that compute the tables in `inst/verification/` are in
`data-raw/benchmark/` in the source repository; its README describes them.

`short_fit.rds` is a fitted object from a short run (five sites,
one chain, 100 draws, `L = 51`, `M = 3`, seed 1234). The tests of the methods
read it, so that they do not need CmdStan. It is not a fit to draw
conclusions from. Its draws are those of the run with version 0.2.0; the
columns `sd` and `map` of its table of site effects were computed again
from those draws with the definitions of version 0.3.0.
