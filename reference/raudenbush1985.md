# Teacher expectancy studies

Results of 19 experiments on whether telling teachers that some of their
pupils are likely to show intellectual growth raises those pupils' IQ
scores. Raudenbush (1984) related the size of the effect to the time
that teachers had known their pupils before the expectation was induced,
and Raudenbush and Bryk (1985) used the studies to illustrate empirical
Bayes meta-analysis. The studies differ widely in size, and so do their
standard errors.

## Usage

``` r
raudenbush1985
```

## Format

A data frame with 19 rows and 11 columns:

- study:

  Study number.

- author, year:

  Authors and year of publication.

- weeks:

  Weeks of contact between the teacher and the pupils before the
  expectation was induced.

- setting:

  Whether the test was given to a group (`"group"`) or individually
  (`"indiv"`).

- tester:

  Whether the person who gave the test knew which pupils had been named
  to the teacher (`"aware"`) or not (`"blind"`).

- n1i, n2i:

  Numbers of pupils in the experimental and the control group.

- yi:

  Standardized mean difference in IQ score. Positive values mean that
  the pupils named to the teacher scored higher on average.

- vi:

  Sampling variance of `yi`.

- sei:

  Standard error of `yi`, the square root of `vi`.

## Source

Raudenbush, S. W. and Bryk, A. S. (1985). Empirical Bayes meta-analysis.
*Journal of Educational Statistics*, 10(2), 75-98.
[doi:10.3102/10769986010002075](https://doi.org/10.3102/10769986010002075)
. The columns `setting` and `tester` are from Raudenbush, S. W. (1984).
Magnitude of teacher expectancy effects on pupil IQ as a function of the
credibility of expectancy induction: A synthesis of findings from 18
experiments. *Journal of Educational Psychology*, 76(1), 85-97.
[doi:10.1037/0022-0663.76.1.85](https://doi.org/10.1037/0022-0663.76.1.85)
. The values were taken from `dat.raudenbush1985` in the metadat
package, which distributes the same data.

## See also

[raudenbush_fit](https://joonho112.github.io/bayesEfron/reference/raudenbush_fit.md),
the model fitted to these data.

## Examples

``` r
head(raudenbush1985)
#>   study             author year weeks setting tester n1i n2i    yi     vi
#> 1     1   Rosenthal et al. 1974     2   group  aware  77 339  0.03 0.0156
#> 2     2        Conn et al. 1968    21   group  aware  60 198  0.12 0.0216
#> 3     3        Jose & Cody 1971    19   group  aware  72  72 -0.14 0.0279
#> 4     4 Pellegrini & Hicks 1972     0   group  aware  11  22  1.18 0.1391
#> 5     5 Pellegrini & Hicks 1972     0   group  blind  11  22  0.26 0.1362
#> 6     6  Evans & Rosenthal 1969     3   group  aware 129 348 -0.06 0.0106
#>         sei
#> 1 0.1249000
#> 2 0.1469694
#> 3 0.1670329
#> 4 0.3729611
#> 5 0.3690528
#> 6 0.1029563
as_bef_data(raudenbush1985)
#> <bef_data>
#> Sites: 19
#> Source: data.frame
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373
```
