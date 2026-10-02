# Extract estimates and standard errors for a fit

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md)
takes a vector of estimates and a vector of their standard errors.
`as_bef_data()` takes the two vectors out of the object that usually
holds them, a list, a data frame or the result of
[`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html),
and checks them. The result is a `bef_data` object, whose elements
`theta_hat` and `sigma` are the arguments of
[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

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
```

## Arguments

- x:

  A list, a data frame, an `escalc` object from the metafor package or a
  `bef_data` object.

- ...:

  For `as_bef_data()`, these dots are for future extensions and must be
  empty. The [`print()`](https://rdrr.io/r/base/print.html) method
  passes them on to [`format()`](https://rdrr.io/r/base/format.html);
  [`format()`](https://rdrr.io/r/base/format.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) do not use them.

- theta_hat, sigma:

  For a data frame, character strings naming the columns of the
  estimates and of the standard errors. `NULL` (the default) means that
  the columns are found by name, as described in Details.

- variance:

  Logical. For a data frame, `TRUE` if the column named in `sigma` holds
  variances; their square roots are then taken. Defaults to `FALSE`.

- site:

  For a data frame, a character string naming a column of site labels.
  Defaults to `NULL`.

- object:

  A `bef_data` object.

- use_cli:

  Logical or `NULL`. With `TRUE` and the cli package installed, the
  first line is set in bold and section titles in color; `FALSE` gives
  plain text. `NULL` (the default) takes the value of
  `getOption("bayesEfron.use_cli", TRUE)`.

## Value

A list of class `bef_data` with the elements `theta_hat`, the estimates;
`sigma`, the standard errors; `names`, the site labels or `NULL`; and
`source`, which is `"list"`, `"data.frame"` or `"metafor::escalc"`
according to the class of `x`.

## Details

The methods differ in where they look for the two vectors.

For a list, the elements `theta_hat` and `sigma` are used. Site labels
are taken from an element `names` or, if there is none, from the names
of `theta_hat`.

For a data frame, the estimates are the column named `theta_hat`, `yi`,
`estimate`, `effect` or `y`, and the standard errors are the column
named `sigma`, `se`, `std_error`, `stderr` or `sei`; in each case the
first of these names that the data frame has. The arguments `theta_hat`
and `sigma` name other columns. If the column holds variances, as the
column `vi` of a metafor data set does, name it in `sigma` and set
`variance = TRUE`. Site labels are taken from the column named in `site`
or, without it, from the row names if they have been set.

For an object of class `escalc`, the estimates are the column `yi` and
the standard errors are the square roots of the column `vi`. Row names
that have been set are used as site labels.

A `bef_data` object is checked again and returned.

At least five sites are needed, and the function stops if a value is
missing or a standard error is not positive. The site labels are shown
when the object is printed; a fitted model numbers the sites in the
order in which they were given.

## Functions

- `print(bef_data)`: Prints the number of sites and the range of the
  estimates and of the standard errors.

- `summary(bef_data)`: Returns a list with the number of sites `K` and
  summaries of `theta_hat` and `sigma`.

- `format(bef_data)`: Returns the lines that
  [`print()`](https://rdrr.io/r/base/print.html) shows, as a character
  vector.

## See also

[`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md),
[raudenbush1985](https://joonho112.github.io/bayesEfron/reference/raudenbush1985.md)

## Examples

``` r
# A data frame: the columns `yi` and `sei` are found by name.
dat <- as_bef_data(raudenbush1985)
dat
#> <bef_data>
#> Sites: 19
#> Source: data.frame
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373

# The same from the variances, with the authors as site labels.
as_bef_data(raudenbush1985, sigma = "vi", variance = TRUE, site = "author")
#> <bef_data>
#> Sites: 19
#> Source: data.frame
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373
#> Names: Rosenthal et al., Conn et al., Jose & Cody, ...

# A list.
as_bef_data(list(theta_hat = raudenbush1985$yi, sigma = raudenbush1985$sei))
#> <bef_data>
#> Sites: 19
#> Source: list
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373

# The two vectors that bayes_efron_fit() takes.
str(dat[c("theta_hat", "sigma")])
#> List of 2
#>  $ theta_hat: num [1:19] 0.03 0.12 -0.14 1.18 0.26 -0.06 -0.02 -0.32 0.27 0.8 ...
#>  $ sigma    : num [1:19] 0.125 0.147 0.167 0.373 0.369 ...

# An escalc object from the metafor package.
es <- metafor::escalc(measure = "GEN", yi = yi, vi = vi, data = raudenbush1985)
as_bef_data(es)
#> <bef_data>
#> Sites: 19
#> Source: metafor::escalc
#> theta_hat: min -0.32; median 0.07; max 1.18
#> sigma: min 0.09381; median 0.167; max 0.373
```
