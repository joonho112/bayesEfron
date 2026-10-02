# Model fitted to the teacher expectancy studies

The object returned by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
for the
[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)
data with the default settings. The examples and the vignettes use it,
so that they run without CmdStan.

## Usage

``` r
raudenbush_fit
```

## Format

An object of class `bef_fit_re`; see
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
for its components.

## Source

Fitted with CmdStan 2.38.0 on 2026-10-02. The script is
`data-raw/raudenbush.R` in the source repository.

## Details

The model was fitted with
`bayes_efron_fit(raudenbush1985$yi, raudenbush1985$sei, seed = 1985)`,
that is, with four chains of 1,000 warmup and 3,000 sampling iterations
each. To keep the package small, one draw in fifteen is stored, so the
object holds 800 of the 12,000 draws, and its posterior summaries are
computed from those 800. The convergence diagnostics that
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
returns are those of the full run. One of the 12,000 iterations reached
the maximum tree depth, which
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
reports as a warning;
[`vignette("diagnostics")`](https://joonho112.github.io/bayesEfron/articles/diagnostics.md)
explains how to read it.

## See also

[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md),
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
raudenbush_fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 19
#> Grid method: paper_realdata
#> Runtime: 1.351 sec
#> Diagnostics: Rhat 1.007; ESS bulk 676.6; ESS tail 888.8; divergences 0; treedepth hits 1
#> Warnings: max_treedepth
#> Use summary() for posterior summaries.
summary(raudenbush_fit)
#> <summary.bef_fit>
#> 
#> Prior g:
#>   mean: 0.1145
#>   var:  0.1094
#>   sd:   0.3009
#> 
#> Diagnostics:
#>   Rhat (max):        1.007
#>   ESS bulk (min):    676.6
#>   ESS tail (min):    888.8
#>   Divergences:       0
#>   Treedepth hits:    1
#>   Effective params:  mean 10.29, sd 2.84
#>   Log marginal lik.: mean -4.047, sd 1.439
#>   Runtime:           1.351 sec
#>   Warnings:          max_treedepth
#> 
#> Theta summary:
#>   site      mean        sd     lower     upper       map
#>   1      0.03805   0.09832     -0.11      0.19      0.04
#>   2       0.0874    0.1121     -0.08      0.28      0.07
#>   3     -0.03442    0.1188     -0.23      0.16     -0.02
#>   4        0.635    0.4979      0.04      1.48      0.19
#>   5       0.1063     0.186     -0.14    0.4015      0.07
#>   6     -0.02282   0.08604     -0.17       0.1     -0.02
#>   14      0.1073    0.1658     -0.11       0.4      0.07
#>   15    -0.05756    0.1177     -0.29    0.1015     -0.05
#>   16    0.001513    0.1169      -0.2      0.19      0.01
#>   17      0.1979    0.1221      0.01       0.4      0.19
#>   18       0.063   0.08101     -0.08      0.19      0.07
#>   19   -0.0006152    0.1197      -0.2      0.19      0.01
#>   ... 7 sites omitted ...
```
