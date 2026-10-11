# Racial callback gaps among 97 large firms

Firm-specific differences in callback rates for applications with
distinctively white and distinctively Black names. The original
correspondence experiment is from Kline, Rose and Walters (2022).
Following their 2024 analysis and the replication code of Walters
(2024), firms with at least 40 sampled vacancies and an overall callback
rate of at least 3 percent are included: 97 firms, 78,910 applications
and 10,453 vacancies. Each firm is one site of the model; the standard
errors are clustered by vacancy within firm.

## Usage

``` r
krw_discrimination
```

## Format

A data frame with 97 rows and 8 columns:

- site:

  Sequential site number, in ascending order of `firm_id`.

- firm_id:

  Firm identifier in the original data.

- firm_name:

  Name of the firm as supplied in the replication data.

- yi:

  Difference in callback probability, white minus Black, estimated by
  ordinary least squares. Positive values favor applications with
  distinctively white names. The values are differences of proportions:
  multiply by 100 for percentage points.

- sei:

  Standard error of `yi`, clustered by vacancy, with the CR1
  finite-sample correction used by Stata's clustered OLS regression.

- vi:

  Sampling variance, `sei^2`.

- n_jobs:

  Number of distinct sampled vacancies in the firm.

- n_obs:

  Number of applications in the firm.

## Source

Kline, P., Rose, E. K. and Walters, C. R. (2022). Systemic
discrimination among large U.S. employers. *The Quarterly Journal of
Economics*, 137(4), 1963-2036.
[doi:10.1093/qje/qjac024](https://doi.org/10.1093/qje/qjac024) . The
analysis sample follows Kline, P., Rose, E. K. and Walters, C. R.
(2024). A discrimination report card. *American Economic Review*,
114(8), 2472-2525.
[doi:10.1257/aer.20230700](https://doi.org/10.1257/aer.20230700) . Data
and regression code were obtained from the replication package
accompanying Walters, C. R. (2024). Empirical Bayes methods in labor
economics. In *Handbook of Labor Economics*, 5, 183-260.
[doi:10.1016/bs.heslab.2024.11.001](https://doi.org/10.1016/bs.heslab.2024.11.001)
. The package is available through
<https://sites.google.com/view/christopher-walters/research>. The
Bayesian reanalysis is described in Lee and Sui (2025),
[doi:10.3390/math13162639](https://doi.org/10.3390/math13162639) .
`data-raw/krw-discrimination.R` in the source repository reproduces the
estimates from `krw_data.dta`, and the attribute `source` of the data
frame records where the data came from.

## See also

[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md),
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)

## Examples

``` r
head(krw_discrimination)
#>   site firm_id          firm_name           yi        sei           vi n_jobs
#> 1    1       1           Rite Aid  0.046957123 0.01619358 0.0002622320    128
#> 2    2       2                KFC  0.022000000 0.01530472 0.0002342345    125
#> 3    3       3  Victoria's Secret  0.042161613 0.02296031 0.0005271759    125
#> 4    4       4     United Rentals  0.005708306 0.01504751 0.0002264274    125
#> 5    5       5    Cardinal Health  0.034077110 0.02145421 0.0004602832    125
#> 6    6       7 Charter / Spectrum -0.010182468 0.01209059 0.0001461824    125
#>   n_obs
#> 1   962
#> 2  1000
#> 3   931
#> 4   917
#> 5   974
#> 6   960
as_bef_data(krw_discrimination)
#> <bef_data>
#> Sites: 97
#> Source: data.frame
#> theta_hat: min -0.02254; median 0.01781; max 0.09817
#> sigma: min 0.005344; median 0.01581; max 0.0392
```
