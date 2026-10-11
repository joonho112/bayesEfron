# Test fixtures

`grid/` holds the inputs and the expected grids for the three grid methods
that follow Lee and Sui (2025). The tests rebuild each grid with
`make_efron_grid()` and compare it to the stored one.

`deconvolveR_baseline.rds` holds empirical Bayes estimates computed with
`deconvolveR::deconv()` for three simulated data sets (K = 20, 50, 100).

`lee_sui_K50.rds` holds the 20 data sets with 50 sites of the benchmark of
Lee and Sui (2025, Appendix E). The files for 100, 200, 500 and 1,500 sites
and the scripts of the benchmark are in `data-raw/benchmark/` in the source
repository; its README describes them.

`grid_comparison_v030.rds` holds comparisons of grids that version 0.3.0
made (`data-raw/grid-comparison.R` in the source repository). A table of
that version has fewer columns and no estimates of g, and the tests read
one to see that it still prints and that `plot()` asks for the comparison
to be run again.

`grid-v030.rds` holds grid points and the values of the spline functions
at them as version 0.3.0 builds them: for 40 simulated sites with each of
the three grid methods of the paper (`"paper_realdata"` with two sets of
arguments), and for the callback data with `"paper_realdata"`. A test
checks that this version builds the same.

`legacy-v01.rds` and `legacy-v02.rds` are fits of the five sites of
`short_fit.rds` in the form in which versions 0.1 and 0.2 stored a fit.
The tests read them to check that such a fit is summarized with the
definitions of this version, with a warning that says so.

`short_fit.rds` is a fitted object from a short run (five sites,
one chain, 100 draws, `L = 51`, `M = 3`, seed 1234). The tests of the methods
read it, so that they do not need CmdStan. It is not a fit to draw
conclusions from. Its draws are those of the run with version 0.2.0; its
table of site effects was computed again from those draws with the
definitions of version 0.4.

`ref-v030.rds` holds draws of version 0.3.0 for two short fits (two
chains, 200 warmup and 100 sampling iterations each): the five sites of
`short_fit.rds` with `L = 51`, `M = 3` and seed 1234, and the 97 firms of
`krw_discrimination` with the default grid and seed 20250807. For each fit
it keeps `alpha`, `lambda`, the summaries of g, `log_g` at every tenth grid
point and the site-level quantities of two sites, with a SHA-256 of every
block of variables of the full draws. It was made from the source
published as 0.3.0, by a script that is not distributed; each fit was run
twice and gave the same draws. The element `conditions` records what the draws
depend on: the Stan file (by its SHA-256), the stanc option `O1`, CmdStan
2.38.0, which writes eight significant digits, the compiler and the
architecture. With any of these different, only agreement within Monte
Carlo error can be expected.
