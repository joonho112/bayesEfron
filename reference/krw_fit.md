# Model fitted to the racial callback gaps

A fitted model for the 97 firms in
[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md),
as returned by
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
with the default settings and with the firm names as site labels. It is
thinned to 800 posterior draws and was created with

    bayes_efron_fit(as_bef_data(krw_discrimination, site = "firm_name"),
                    seed = 20250807, thin = 15)

## Usage

``` r
krw_fit
```

## Format

An object of class `bef_fit_re`; see
[bef_fit](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
for its components and methods.

## Source

Fitted with CmdStan 2.38.0 on 2026-10-10. The script is
`data-raw/example-fits.R` in the source repository. The sources of the
data are in
[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md).

## Details

Each of the 4 chains ran 1,000 warmup and 3,000 sampling iterations, and
one draw out of every 15 was stored. Summaries are computed from the 800
stored draws, whereas the sampler diagnostics reported by
[`diagnose()`](https://joonho112.github.io/bayesEfron/reference/diagnose.md)
are computed from all 12,000. For these data the `"pilot"` method does
not narrow the range: the grid has 101 points and extends beyond the
range of the estimates by half the width of that range on each side. The
object contains the draws of \\g\\ but no site-level draws, and it can
be used without CmdStan.

Sites are in the same order as the rows of
[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md).
Tables, plots and the names of the
[`coef()`](https://rdrr.io/r/stats/coef.html) vector use the firm names
as site labels, and a site can be selected by position or by name.

## See also

[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md),
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)

## Examples

``` r
krw_fit
#> <bayesEfron fit>
#> Model family: RE
#> Sites: 97
#> Grid method: pilot
#> Runtime: 0.7197 sec
#> Draws: 800 (4 chains, one iteration in 15 kept); site-level draws not kept
#> Diagnostics: Rhat 1.002; ESS bulk 1731; ESS tail 1780; divergences 0; treedepth hits 0
#> Standard errors: 0.005344 to 0.0392 (largest / smallest 7.335)
#> Grid: 101 points from -0.0829 to 0.1585; smallest standard error 2.213 grid spacings
#> Grid range: rule "pilot", not narrowed (the range of "paper_realdata"); point masses; 6 spline functions
#> Next: summary() and site_summary() for the sites, g_summary() for the distribution of effects, plot().
confint(krw_fit, type = "g")
#> 90% posterior intervals (equal-tailed)
#>   parameter         mean        lower        upper   site        point
#> 1    mean_g 0.0166017761 0.0130811925 0.0206491286 mean_g 0.0166017761
#> 2     var_g 0.0003076948 0.0001399332 0.0005738865  var_g 0.0003076948
#> 3      sd_g 0.0171621568 0.0118293270 0.0239559271   sd_g 0.0171621568
site_summary(krw_fit, sites = c("KFC", "Rite Aid"), scale = 100)
#>   site    label estimate       se     mean       sd     lower    upper     mode
#> 1    2      KFC 2.200000 1.530472 1.786494 1.052594 0.1227187 3.581179 1.608692
#> 2    1 Rite Aid 4.695712 1.619358 3.048175 1.296795 1.0775168 5.349968 2.815804
#>   sd_given_g sd_about_g
#> 1   1.040999  0.1558017
#> 2   1.262179  0.2976259
```
