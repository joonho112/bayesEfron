# Extract estimates and standard errors from a data frame or list

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
needs a vector of effect estimates and a vector of their standard
errors. `as_bef_data()` extracts the two vectors from a data frame, a
list or an `escalc` object from the metafor package and checks them. The
result is a `bef_data` object, which can be passed to
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
in place of the two vectors.

## Usage

``` r
as_bef_data(x, ...)

# Default S3 method
as_bef_data(x, ...)

# S3 method for class 'list'
as_bef_data(x, ...)

# S3 method for class 'escalc'
as_bef_data(x, ...)

# S3 method for class 'data.frame'
as_bef_data(
  x,
  ...,
  theta_hat = NULL,
  sigma = NULL,
  variance = FALSE,
  site = NULL
)

# S3 method for class 'bef_data'
print(x, ...)

# S3 method for class 'bef_data'
summary(object, ...)

# S3 method for class 'bef_data'
format(x, ..., use_cli = NULL)

# S3 method for class 'bef_data'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  For `as_bef_data()`, a data frame, a list, an `escalc` object from the
  metafor package or a `bef_data` object; see Details. For
  [`print()`](https://rdrr.io/r/base/print.html),
  [`format()`](https://rdrr.io/r/base/format.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), a
  `bef_data` object.

- ...:

  For `as_bef_data()`, these dots are for future extensions and must be
  empty. The [`print()`](https://rdrr.io/r/base/print.html) method
  passes them on to [`format()`](https://rdrr.io/r/base/format.html);
  the other methods ignore them.

- theta_hat, sigma:

  For a data frame, the names of the columns that contain the estimates
  and the standard errors, as character strings. If `NULL` (the
  default), the columns are found by name; see Details.

- variance:

  For a data frame, set to `TRUE` if the column specified in `sigma`
  contains variances. Their square roots are then used as the standard
  errors.

- site:

  For a data frame, an optional character string naming a column of site
  labels. If `NULL` (the default), row names that have been set are used
  as labels.

- object:

  A `bef_data` object.

- use_cli:

  For [`format()`](https://rdrr.io/r/base/format.html): if `TRUE` and
  the cli package is installed, the first line is in bold; if `FALSE`,
  the result is plain text. If `NULL` (the default), the value of
  `getOption("bayesEfron.use_cli", TRUE)` is used.

- row.names, optional:

  Not used; they are arguments of the generic function.

## Value

An object of class `bef_data`, which is a list with the following
components.

- `theta_hat`:

  The estimates.

- `sigma`:

  The standard errors.

- `names`:

  The site labels, or `NULL` if there are none.

- `source`:

  `"list"`, `"data.frame"` or `"metafor::escalc"`, depending on the
  class of `x`.

## Details

The methods differ in where they look for the estimates and the standard
errors.

For a list, the elements `theta_hat` and `sigma` are used. Site labels
are taken from an element `names` or, if there is none, from the names
of `theta_hat`.

For a data frame, the estimates are the column named `theta_hat`, `yi`,
`estimate`, `effect` or `y`, and the standard errors are the column
named `sigma`, `se`, `std_error`, `stderr` or `sei`. If the data frame
has more than one of these columns, the first in the list is used. If
the data frame has none of the columns for standard errors but a column
`vi`, `var` or `variance`, as a metafor data set does, the square roots
of that column are used, with a message. Use the arguments `theta_hat`
and `sigma` to specify columns with other names. For a column of
variances with another name, specify it as `sigma` and set
`variance = TRUE`. If `sigma` names a column whose name begins with `vi`
or `var` and `variance` is not `TRUE`, a warning is issued, because the
name suggests that the column contains variances. Site labels are taken
from the column named in `site`. Without `site`, the row names are used
if they have been set.

For an `escalc` object, the estimates are the column `yi` and the
standard errors are the square roots of the column `vi`. Site labels are
the study labels given to
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
in its argument `slab` or, without them, row names that have been set.

A `bef_data` object is checked again and returned.

At least 5 sites are required. An error is thrown if a value is missing,
if a standard error is not positive, or if a site label occurs more than
once. The site labels are stored in a fitted model: they are shown in
its tables and plots, and sites can be selected by their labels.

## Functions

- `print(bef_data)`: Prints the number of sites and the ranges of the
  estimates and of the standard errors.

- `summary(bef_data)`: Returns a list with the number of sites `K` and
  summaries of `theta_hat` and `sigma`.

- `format(bef_data)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows as a character
  vector.

- `as.data.frame(bef_data)`: Returns a data frame with one row per site
  and the columns `site`, `label`, `estimate` and `se`. If there are no
  site labels, `label` is the site number.

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[krw_discrimination](https://joonho112.github.io/bayesEfron/reference/krw_discrimination.md)

## Examples

``` r
# A data frame: the columns `yi` and `sei` are found by name.
dat <- as_bef_data(krw_discrimination)
dat
#> <bef_data>
#> Sites: 97
#> Source: data.frame
#> theta_hat: min -0.02254; median 0.01781; max 0.09817
#> sigma: min 0.005344; median 0.01581; max 0.0392

# The same from the variances, with firm IDs as site labels.
as_bef_data(krw_discrimination, sigma = "vi", variance = TRUE, site = "firm_id")
#> <bef_data>
#> Sites: 97
#> Source: data.frame
#> theta_hat: min -0.02254; median 0.01781; max 0.09817
#> sigma: min 0.005344; median 0.01581; max 0.0392
#> Names: 1, 2, 3, ...

# A list.
as_bef_data(list(theta_hat = krw_discrimination$yi, sigma = krw_discrimination$sei))
#> <bef_data>
#> Sites: 97
#> Source: list
#> theta_hat: min -0.02254; median 0.01781; max 0.09817
#> sigma: min 0.005344; median 0.01581; max 0.0392

# The two vectors that bayes_efron_fit() takes.
str(dat[c("theta_hat", "sigma")])
#> List of 2
#>  $ theta_hat: num [1:97] 0.04696 0.022 0.04216 0.00571 0.03408 ...
#>  $ sigma    : num [1:97] 0.0162 0.0153 0.023 0.015 0.0215 ...

# An escalc object from the metafor package.
es <- metafor::escalc(measure = "GEN", yi = yi, vi = vi, data = krw_discrimination)
as_bef_data(es)
#> <bef_data>
#> Sites: 97
#> Source: metafor::escalc
#> theta_hat: min -0.02254; median 0.01781; max 0.09817
#> sigma: min 0.005344; median 0.01581; max 0.0392
```
