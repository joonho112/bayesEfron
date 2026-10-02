# Using metafor and data frames as input

The model needs two numbers for each site: an estimate of the effect and
the standard error of that estimate. Data seldom arrive in that form. A
table may hold variances where standard errors are wanted, or the
estimates may still have to be computed from the means, standard
deviations and sample sizes of two groups. The sites usually have names
as well, which you will want to see again next to the results.
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
deals with the first two points and keeps the names.

## Estimates already in a data frame

The teacher expectancy studies of Raudenbush and Bryk
([1985](#ref-raudenbush1985)) come as a data frame with one row for each
study.

``` r

studies <- raudenbush1985[, c("author", "yi", "vi", "sei")]
studies$sei <- round(studies$sei, 4)
head(studies, 4)
#>               author    yi     vi    sei
#> 1   Rosenthal et al.  0.03 0.0156 0.1249
#> 2        Conn et al.  0.12 0.0216 0.1470
#> 3        Jose & Cody -0.14 0.0279 0.1670
#> 4 Pellegrini & Hicks  1.18 0.1391 0.3730
```

`yi` is the estimate, a standardized mean difference, `vi` is its
sampling variance and `sei` its standard error.
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
finds the two columns it needs by their names.

``` r

dat <- as_bef_data(raudenbush1985)
dat
#> <bef_data>
#> Sites: 19
#> Source: data.frame
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373
```

For the estimates it takes the first of the names `theta_hat`, `yi`,
`estimate`, `effect` and `y` that the data frame has, and for the
standard errors the first of `sigma`, `se`, `std_error`, `stderr` and
`sei`. Columns with other names are given in the arguments `theta_hat`
and `sigma`. If your second column holds variances, as `vi` does, set
`variance = TRUE` and the function takes the square roots.

``` r

dat_v <- as_bef_data(raudenbush1985, theta_hat = "yi", sigma = "vi",
                     variance = TRUE)
all.equal(dat_v$sigma, dat$sigma)
#> [1] TRUE
```

The function stops if a value is missing, if a standard error is not
positive or if there are fewer than five sites. It does not drop
incomplete rows, because the number of sites would then change without
notice.

``` r

incomplete <- raudenbush1985
incomplete$sei[3] <- NA
try(as_bef_data(incomplete))
#> Error in as_bef_data(incomplete) : 
#>   Columns `yi` and `sei` must not contain missing values (found 0 and 1).
```

## From group summaries to estimates

Often the data are one step further back, and each study reports the
mean, the standard deviation and the size of two groups.
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
computes an estimate and its sampling variance from such summaries
([Viechtbauer 2010](#ref-viechtbauer2010)), and
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
accepts its result. The example is the nine studies of the length of
hospital stay of stroke patients under specialized and under routine
care that Normand ([1999](#ref-normand1999)) analyzed. The data are in
the metadat package, which metafor loads.

``` r

library(metafor)
stroke <- escalc(
  measure = "MD",
  m1i = m1i, sd1i = sd1i, n1i = n1i,
  m2i = m2i, sd2i = sd2i, n2i = n2i,
  data = dat.normand1999
)
stroke[, c("source", "m1i", "m2i", "yi", "vi")]
#> 
#>               source m1i m2i       yi       vi 
#> 1          Edinburgh  55  75 -20.0000  40.5080 
#> 2     Orpington-Mild  27  29  -2.0000   2.0806 
#> 3 Orpington-Moderate  64 119 -55.0000  15.6984 
#> 4   Orpington-Severe  66 137 -71.0000 150.2222 
#> 5      Montreal-Home  14  18  -4.0000  17.3077 
#> 6  Montreal-Transfer  19  18   1.0000   1.1673 
#> 7          Newcastle  52  41  11.0000  94.5891 
#> 8               Umea  21  31 -10.0000   6.3109 
#> 9            Uppsala  30  23   7.0000  19.8423
```

With `measure = "MD"` the estimate `yi` is the difference between the
two means, here in days, and `vi` is the sum of the two squared standard
errors of the means, \\s_1^2 / n_1 + s_2^2 / n_2\\. A difference in the
original units is the natural choice when every study measures the
outcome in the same way, as these do. When the studies use different
instruments, as the teacher expectancy studies did with their IQ tests,
the standardized mean difference (`measure = "SMD"`) puts the estimates
on a common scale.

``` r

rownames(stroke) <- stroke$source
stroke_dat <- as_bef_data(stroke)
stroke_dat
#> <bef_data>
#> Sites: 9
#> Source: metafor::escalc
#> theta_hat: min -71; median -4; max 11
#> sigma: min 1.08; median 4.16; max 12.26
#> Names: Edinburgh, Orpington-Mild, Orpington-Moderate, ...
```

For an `escalc` object the estimates are the column `yi` and the
standard errors are the square roots of `vi`. The site labels are taken
from the row names, which is why the row names were set first; labels
given to
[`escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
in its argument `slab` are not used.

Nine studies are enough for the function to run, but they can say little
about the shape of the distribution of effects, and with so few the
results depend more on the choices that
[`vignette("choosing-a-grid")`](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md)
describes. The check of the intervals that
[`vignette("implementation")`](https://joonho112.github.io/bayesEfron/articles/implementation.md)
reports starts at 50 sites; we have no evidence on how they perform with
fewer.

## A list

When the estimates and standard errors are two vectors, they can be
passed to
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
directly.
[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md)
accepts them as a list as well, which is a way to have them checked
beforehand.

``` r

as_bef_data(list(
  theta_hat = c(0.03, 0.12, -0.14, 1.18, 0.26),
  sigma = c(0.125, 0.147, 0.167, 0.373, 0.369),
  names = c("Rosenthal", "Conn", "Jose", "Pellegrini 1", "Pellegrini 2")
))
#> <bef_data>
#> Sites: 5
#> Source: list
#> theta_hat: min -0.14; median 0.12; max 1.18
#> sigma: min 0.125; median 0.167; max 0.373
#> Names: Rosenthal, Conn, Jose, ...
```

## Fitting, and putting the names back

The elements `theta_hat` and `sigma` of the result are the two arguments
of
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

``` r

dat <- as_bef_data(raudenbush1985, site = "author")
fit <- bayes_efron_fit(dat$theta_hat, dat$sigma, seed = 1985)
```

The fitted object numbers the sites in the order of the data and does
not keep their labels, so the labels have to be put back by hand. The
object `raudenbush_fit` in the package holds the result of this fit and
keeps one posterior draw in fifteen.

``` r

dat <- as_bef_data(raudenbush1985, site = "author")
fit <- raudenbush_fit
results <- cbind(study = dat$names, as.data.frame(fit))
results$mean <- round(results$mean, 2)
head(results[, c("study", "mean", "hpdi_lower", "hpdi_upper")])
#>                study  mean hpdi_lower hpdi_upper
#> 1   Rosenthal et al.  0.04      -0.11     0.1900
#> 2        Conn et al.  0.09      -0.08     0.2800
#> 3        Jose & Cody -0.03      -0.23     0.1600
#> 4 Pellegrini & Hicks  0.63       0.04     1.4800
#> 5 Pellegrini & Hicks  0.11      -0.14     0.4015
#> 6  Evans & Rosenthal -0.02      -0.17     0.1000
```

`mean` is the posterior mean of the effect. `hpdi_lower` and
`hpdi_upper` are the 5% and 95% posterior quantiles, the limits of a 90%
interval; in spite of their names the interval is equal-tailed.
`confint(fit, level = 0.95)` gives 95% intervals, the level that metafor
uses by default for its confidence intervals, and
[`?bef_fit`](https://joonho112.github.io/bayesEfron/reference/bef_fit.md)
describes all the columns.

## What differs from a random-effects model

metafor fits the random-effects model to the same two columns. The
difference between the two analyses lies in what is assumed about the
distribution of true effects. The random-effects model takes it to be
normal, with a mean and a standard deviation that are estimated.
bayesEfron estimates the distribution itself, \\g\\, and reports its
mean and standard deviation as summaries.

``` r

normal_fit <- rma(yi, vi, data = raudenbush1985)
c(mean = unname(coef(normal_fit)), sd = sqrt(normal_fit$tau2))
#>       mean         sd 
#> 0.08370824 0.13720914
confint(fit, type = "g")[c(1, 3), ]
#>     site       lower     upper     point
#> 1 mean_g 0.003896369 0.2499901 0.1144608
#> 3   sd_g 0.104426290 0.5515191 0.3008526
```

The two means are close, whereas the standard deviation of \\g\\ is
about twice that of the normal model. Both numbers describe how much the
true effects vary between studies, and they differ because the two
models allow different shapes for that variation. Part of the difference
is the long right tail of the estimated \\g\\ (it is shown in
[`vignette("bayesEfron")`](https://joonho112.github.io/bayesEfron/articles/bayesEfron.md)),
which a normal distribution cannot have. The other part is that \\g\\
keeps a little probability beyond the smallest and the largest estimate,
where the data hardly bear on it: 0.06 in all, which accounts for 59% of
the variance of the estimated distribution, the posterior mean of \\g\\.
A standard deviation is sensitive to such probability, so that of \\g\\
says as much about the tails that the model allows as about the bulk of
the effects.

The site effects differ too.
[`blup()`](https://wviechtb.github.io/metafor/reference/blup.html) gives
the prediction of the normal model for each study, which is drawn toward
the overall mean by a factor that depends only on the standard error.
With bayesEfron, how far an estimate is drawn in depends on the shape of
\\g\\ as well.

``` r

compare <- data.frame(
  study = raudenbush1985$author,
  estimate = raudenbush1985$yi,
  normal = round(blup(normal_fit)$pred, 2),
  bayesEfron = round(coef(fit), 2)
)
compare[order(-compare$estimate)[1:6], ]
#>                   study estimate normal bayesEfron
#> 4    Pellegrini & Hicks     1.18   0.21       0.63
#> 10              Maxwell     0.80   0.25       0.47
#> 11               Carter     0.54   0.16       0.22
#> 17 Rosenthal & Jacobson     0.30   0.19       0.20
#> 9                Kester     0.27   0.16       0.16
#> 5    Pellegrini & Hicks     0.26   0.11       0.11
```

For most of the studies the two columns are close, but they part for the
largest estimates. Under the normal model the largest estimate, which
has a large standard error, is taken to be mostly noise and is drawn
most of the way to the mean. Under the estimated \\g\\, with its long
right tail, a true effect of that size is less surprising, and the
estimate is drawn in less far. Which of the two is nearer the truth
cannot be read from 19 studies. The comparison shows that the answer for
the extreme studies depends on what is assumed about the distribution.

## References

Normand, Sharon-Lise T. 1999. “Meta-Analysis: Formulating, Evaluating,
Combining, and Reporting.” *Statistics in Medicine* 18 (3): 321–59.
<https://doi.org/10.1002/(SICI)1097-0258(19990215)18:3%3C321::AID-SIM28%3E3.0.CO;2-P>.

Raudenbush, Stephen W., and Anthony S. Bryk. 1985. “Empirical Bayes
Meta-Analysis.” *Journal of Educational Statistics* 10 (2): 75–98.
<https://doi.org/10.3102/10769986010002075>.

Viechtbauer, Wolfgang. 2010. “Conducting Meta-Analyses in R with the
metafor Package.” *Journal of Statistical Software* 36 (3).
<https://doi.org/10.18637/jss.v036.i03>.
