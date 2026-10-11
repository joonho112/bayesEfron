# Small-class effects in 79 schools of Project STAR

School-specific estimates of the effect of small classes on achievement
in kindergarten, from Tennessee's Student Teacher Achievement Ratio
(STAR) project. Within each of 79 schools, students and teachers were
assigned at random to small classes or to regular classes, the regular
ones with or without a teacher's aide, so the experiment is a multisite
trial with the schools as its sites. The estimate of a school is the
difference between its small classes and its regular classes of either
kind on a composite of the kindergarten reading and mathematics scores.

## Usage

``` r
star_schools
```

## Format

A data frame with 79 rows and 12 columns:

- site:

  Sequential site number, in ascending order of `school_id`.

- school_id:

  School identifier in the public-use data.

- location:

  Kind of area that the school serves: `"inner city"`, `"suburban"`,
  `"rural"` or `"urban"`.

- yi:

  Estimate of the effect of small classes: small minus regular, in
  standard deviations of the composite score.

- sei:

  Standard error of `yi`, from the model with a random intercept for
  each class.

- vi:

  Sampling variance, `sei^2`.

- yi_ols, sei_ols:

  The same contrast and its standard error by ordinary least squares
  within the school.

- n_students, n_students_small:

  Numbers of students with both scores: in the school, and in its small
  classes.

- n_classes, n_classes_small:

  Numbers of classes: in the school, and small ones.

## Source

Achilles, C. M., Bain, H. P., Bellott, F., Boyd-Zaharias, J., Finn, J.,
Folger, J., Johnston, J. and Word, E. (2008). Tennessee's Student
Teacher Achievement Ratio (STAR) project. Harvard Dataverse, V1.
[doi:10.7910/DVN/SIWH9F](https://doi.org/10.7910/DVN/SIWH9F) . The data
are released under CC0 1.0. The experiment is described in Finn, J. D.
and Achilles, C. M. (1990). Answers and questions about class size: A
statewide experiment. *American Educational Research Journal*, 27(3),
557-577.
[doi:10.3102/00028312027003557](https://doi.org/10.3102/00028312027003557)
. `data-raw/star.R` in the source repository builds the estimates from
the student file `STAR_Students.tab`, and the attribute `source` of the
data frame records the file, the sample and the model.

## Details

The outcome is the mean of the standardized reading and mathematics
scale scores of the Stanford Achievement Test, standardized again, so
the estimates are in standard deviations of that composite. Of the 6,325
students in kindergarten, the 5,786 with both scores are used; they sat
in 325 classes, 127 of them small.

The estimates and their standard errors come from one linear mixed model
fitted to all the schools by restricted maximum likelihood. The model
has an intercept and a small-class contrast for each school and a random
intercept for each class. The variance of the class intercepts is taken
to be the same for small and for regular classes and in every school,
and so is the residual variance, so the standard error of a school
depends only on the numbers and the sizes of its small and its regular
classes. A school has from one to four small classes, and that number,
more than the number of its students, limits how precise its estimate
can be.

`yi_ols` and `sei_ols` are the same contrast computed within each school
by ordinary least squares, with the students of a school taken as
independent. That ignores what the students of a class have in common,
and the standard error is smaller for every school. The two columns are
there to show how much a fit depends on its standard errors; the article
"Small classes in 79 schools" on the package website,
<https://joonho112.github.io/bayesEfron/>, fits both, and it also varies
the assumption about the variance of the classes.

## See also

[`as_bef_data()`](https://joonho112.github.io/bayesEfron/reference/as_bef_data.md),
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)

## Examples

``` r
head(star_schools)
#>   site school_id location          yi       sei        vi      yi_ols   sei_ols
#> 1    1    112038    urban -0.22151425 0.4924872 0.2425436 -0.22519182 0.1344110
#> 2    2    123056    rural  0.07580528 0.4930496 0.2430979  0.07736164 0.1712559
#> 3    3    128068    rural -0.02022031 0.4980018 0.2480058 -0.02650660 0.2273065
#> 4    4    128076 suburban  0.05546138 0.4962996 0.2463133  0.05676752 0.1630621
#> 5    5    128079    rural  1.39287138 0.4894715 0.2395823  1.40446846 0.2903351
#> 6    6    130085 suburban  0.15082600 0.3991840 0.1593479  0.14366035 0.1824108
#>   n_students n_students_small n_classes n_classes_small
#> 1         52               15         3               1
#> 2         51               15         3               1
#> 3         54               13         3               1
#> 4         58               13         3               1
#> 5         59               15         3               1
#> 6         71               29         4               2
as_bef_data(star_schools)
#> <bef_data>
#> Sites: 79
#> Source: data.frame
#> theta_hat: min -0.8986; median 0.1459; max 1.393
#> sigma: min 0.2843; median 0.4081; max 0.5684

# The contrast and standard errors that take the students of a school as
# independent.
as_bef_data(star_schools, theta_hat = "yi_ols", sigma = "sei_ols")
#> <bef_data>
#> Sites: 79
#> Source: data.frame
#> theta_hat: min -0.8664; median 0.1437; max 1.404
#> sigma: min 0.1216; median 0.2237; max 0.3945
```
