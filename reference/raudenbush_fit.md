# Model fitted to the teacher expectancy studies

A fitted model for the 19 studies in
[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md),
as returned by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
with the default settings. It is thinned to 800 posterior draws, is used
in
[`vignette("teacher-expectancy")`](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md),
and was created with

    bayes_efron_fit(raudenbush1985, seed = 1985, thin = 15)

## Usage

``` r
raudenbush_fit
```

## Format

An object of class `bef_fit_re`; see
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
for its components and methods.

## Source

Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
`data-raw/example-fits.R` in the source repository.

## Details

Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations, and
one draw out of every 15 was stored. Summaries are computed from the 800
stored draws, whereas the sampler diagnostics reported by
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
are computed from all 12,000. For these 19 studies the `"pilot"` method
does not narrow the range, so the grid is the same as with
`grid_method = "paper_realdata"`. The object contains the draws of \\g\\
but no site-level draws, and it can be used without CmdStan.

## See also

[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md),
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[`vignette("teacher-expectancy")`](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md)

## Examples

``` r
raudenbush_fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 19
#> Grid method: pilot
#> Runtime: 0.6065 sec
#> Draws: 800 (4 chains, one iteration in 15 kept); site-level draws not kept
#> Diagnostics: Rhat 1.003; ESS bulk 994.5; ESS tail 1185; divergences 0; treedepth hits 0
#> Standard errors: 0.09381 to 0.373 (largest / smallest 3.976)
#> Grid: 101 points from -1.07 to 1.93; smallest standard error 3.127 grid spacings
#> Grid range: rule "pilot", not narrowed (the range of "paper_realdata"); point masses; 6 spline functions
#> Next: summary() and site_summary() for the sites, g_summary() for the distribution of effects, plot().
summary(raudenbush_fit)
#> <summary.bef_fit>
#> 
#> Prior g:
#>   mean: 0.1167
#>   var:  0.1159
#>   sd:   0.3123
#>   (each is a posterior mean, so sd squared is not var)
#> 
#> Diagnostics:
#>   Rhat (max):        1.003
#>   ESS bulk (min):    994.5
#>   ESS tail (min):    1185
#>   Divergences:       0
#>   Treedepth hits:    0
#>   Effective params:  mean 10.62, sd 2.655
#>   Log marginal lik.: mean -4.168, sd 1.524
#>   Runtime:           0.6065 sec
#> 
#> Site effects (posterior summaries, 90% intervals):
#>   site      mean        sd     lower     upper      mode
#>   1      0.03778    0.0993   -0.1227    0.2033      0.04
#>   2       0.0883    0.1134  -0.09111    0.2822      0.07
#>   3      -0.0366    0.1198   -0.2388    0.1556     -0.02
#>   4       0.6578    0.5102   0.02838     1.611      0.19
#>   5       0.1079    0.1881   -0.1574    0.4363      0.07
#>   6     -0.02424   0.08638   -0.1686    0.1163     -0.02
#>   ... 7 sites omitted ...
#>   14      0.1091    0.1681   -0.1377    0.4036      0.07
#>   15    -0.06035    0.1184   -0.2631     0.127     -0.05
#>   16   0.0003407    0.1182   -0.1948    0.1949      0.01
#>   17      0.2009    0.1229  0.008355     0.412      0.19
#>   18     0.06311   0.08166  -0.06939    0.1999      0.07
#>   19    -0.00186     0.121   -0.2022    0.1968      0.01
```
