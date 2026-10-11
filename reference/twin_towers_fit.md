# Model fitted to the Twin Towers simulated data

A fitted model for the 1,500 simulated sites in
[twin_towers](https://joonho112.github.io/bayesEfron/reference/twin_towers.md),
as returned by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
with 8 spline functions. The true effects, which are also in the data
frame, are not used in the fit. The fit is thinned to 800 posterior
draws and was created with

    bayes_efron_fit(twin_towers, M = 8, seed = 20250627, thin = 15)

## Usage

``` r
twin_towers_fit
```

## Format

An object of class `bef_fit_re`; see
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
for its components and methods.

## Source

Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
`data-raw/example-fits.R` in the source repository. The sources of the
data are in
[twin_towers](https://joonho112.github.io/bayesEfron/reference/twin_towers.md).

## Details

In this data set, 8 spline functions (`M = 8`) give an estimate of \\g\\
that is closer to the true effects than the default of 6, and, unlike
with 6, no iteration reaches the maximum tree depth. See
[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)
for a comparison of the two settings, including the lower coverage of
the site intervals with 8 functions in simulated data sets. For these
data the `"pilot"` method narrows the range of `"paper_realdata"`, and
the printed fit shows by how much.

Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations, and
one draw out of every 15 was stored. Summaries are computed from the 800
stored draws, whereas the sampler diagnostics reported by
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
are computed from all 12,000. The object contains the draws of \\g\\ but
no site-level draws, and it can be used without CmdStan. Sites are
numbered in the order of the data.

A comparison of this fit with fits of the same data on other ranges is
stored in the file `twin_towers_comparison.rds` in the `examples`
directory of the installed package. The other ranges include the range
around the true effects (`grid_method = "paper_simulation"`), which is
used in the simulations of Lee and Sui (2025), and the wide range of
`"paper_realdata"`.

## See also

[twin_towers](https://joonho112.github.io/bayesEfron/reference/twin_towers.md),
[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)

## Examples

``` r
twin_towers_fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 1500
#> Grid method: pilot
#> Runtime: 12.41 sec
#> Draws: 800 (4 chains, one iteration in 15 kept); site-level draws not kept
#> Diagnostics: Rhat 1.002; ESS bulk 1287; ESS tail 993.3; divergences 0; treedepth hits 0
#> Standard errors: 0.2182 to 1.964 (largest / smallest 9)
#> Grid: 101 points from -3.003 to 3.674; smallest standard error 3.268 grid spacings
#> Grid range: rule "pilot", narrowed to 30% of the range of "paper_realdata"; point masses; 8 spline functions
#> Next: summary() and site_summary() for the sites, g_summary() for the distribution of effects, plot().
confint(twin_towers_fit, type = "g")
#> 90% posterior intervals (equal-tailed)
#>   parameter      mean     lower     upper   site     point
#> 1    mean_g 0.7371258 0.6632461 0.8090581 mean_g 0.7371258
#> 2     var_g 2.1731465 2.0598402 2.2752720  var_g 2.1731465
#> 3      sd_g 1.4739846 1.4352143 1.5084005   sd_g 1.4739846
```
