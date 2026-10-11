# Twin Towers simulated data

A simulated data set of 1,500 sites from the public replication data of
Lee and Sui (2025). The true effects are the 1,500 values of
`disjointTheta` in the deconvolveR package, which lie in two separated
groups; normally distributed sampling errors obscure the two groups in
the estimates. The largest standard error is nine times the smallest
(`R = 9`), and the design parameter `I` is 0.7; both are defined below.

## Usage

``` r
twin_towers
```

## Format

A data frame with 1,500 rows and 7 columns:

- site:

  Site number, in the order of the replication data.

- theta_true:

  True effect, a value of `disjointTheta`.

- yi:

  Estimate of the effect.

- sei:

  Standard error of the estimate.

- vi:

  Sampling variance, `sei^2`.

- I:

  Design parameter, \\1 / (1 + \bar v)\\ with \\\bar v\\ the geometric
  mean of `vi`. It would be the reliability of the estimates if the true
  effects had variance one; their variance in this data set is 2.11.

- R:

  Ratio of the largest to the smallest standard error. The ratio of the
  variances is `R^2`, 81 in this data set.

## Source

Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
2639. [doi:10.3390/math13162639](https://doi.org/10.3390/math13162639) .
Public replication data:
<https://github.com/joonho112/fully-bayes-efron-prior>. The true effects
are from the Twin Towers example distributed in deconvolveR; see
Narasimhan, B. and Efron, B. (2020). deconvolveR: A G-modeling program
for deconvolution and empirical Bayes estimation. *Journal of
Statistical Software*, 94(11), 1-20.
[doi:10.18637/jss.v094.i11](https://doi.org/10.18637/jss.v094.i11) .
`data-raw/twin-towers.R` in the source repository selects the data set
and checks it, and the attribute `source` of the data frame records
where it came from.

## Details

The sampling variances are equally spaced on a log scale and were
assigned to the sites at random before the normal errors were added. The
data set is taken as it stands from the replication data of that study,
where it is called a panel.

## See also

[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md),
[`vignette("twin-towers")`](https://joonho112.github.io/bayesEfron/articles/twin-towers.md)

## Examples

``` r
head(twin_towers)
#>   site theta_true        yi       sei         vi   I R
#> 1    1   1.575405 1.8642862 1.7828639 3.17860355 0.7 9
#> 2    2   1.474004 1.4128779 1.1268516 1.26979464 0.7 9
#> 3    3   1.483298 1.1633946 0.2613307 0.06829374 0.7 9
#> 4    4   1.058912 1.7657197 1.4865568 2.20985114 0.7 9
#> 5    5   1.315141 0.8854123 0.3304032 0.10916626 0.7 9
#> 6    6   1.222705 0.9433596 0.2857737 0.08166658 0.7 9
as_bef_data(twin_towers)
#> <bef_data>
#> Sites: 1500
#> Source: data.frame
#> theta_hat: min -5.118; median 1.067; max 6.01
#> sigma: min 0.2182; median 0.6547; max 1.964
```
