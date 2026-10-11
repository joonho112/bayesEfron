# Export the estimated distribution of effects for the ebnm package

`as_normalmix()` returns the posterior mean of \\g\\ as a `normalmix`
object, the class that the ashr and ebnm packages use for a prior
distribution that is a mixture of normal distributions. The means of the
mixture components are the grid points. The ebnm package can use this
object to compute the posterior distributions of the effects with \\g\\
fixed at that estimate, so the fit can be compared with empirical Bayes
fits of the same data.

## Usage

``` r
as_normalmix(x, ...)
```

## Arguments

- x:

  A `bef_fit` object created by
  [`bayes_efron_fit()`](https://joonho112.github.io/bayesEfron/reference/bayes_efron_fit.md).

- ...:

  These dots are for future extensions and must be empty.

## Value

An object of class `normalmix`, as created by
[`ashr::normalmix()`](https://rdrr.io/pkg/ashr/man/normalmix.html). It
is a list with the following components, each with one value per grid
point.

- `pi`:

  The mixture proportions.

- `mean`:

  The component means, which are the grid points.

- `sd`:

  The component standard deviations.

## Details

The mixture has one component per grid point. The mixture proportions
are the posterior means of the probabilities of \\g\\. The component
standard deviations are zero for a fit with point masses. For a fit with
`component = "normal"`, they all equal the standard deviation of its
normal components.

Pass the object to
[`ebnm::ebnm_npmle()`](https://rdrr.io/pkg/ebnm/man/ebnm_npmle.html) as
`g_init` with `fix_g = TRUE` to obtain the plug-in posterior
distributions of the effects, in which \\g\\ is treated as known. They
differ from the posterior distributions of the fit, which average over
the posterior draws of \\g\\. For most sites the posterior means are
close and the standard deviation from the fit is slightly larger, but
neither is true for every site. See the Details of
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md)
for the relation between the two standard deviations.

`as_normalmix()` constructs the object itself, so neither ashr nor ebnm
needs to be installed.

## See also

[`g_summary()`](https://joonho112.github.io/bayesEfron/reference/g_summary.md),
[`site_summary()`](https://joonho112.github.io/bayesEfron/reference/site_summary.md),
[`vignette("empirical-bayes")`](https://joonho112.github.io/bayesEfron/articles/empirical-bayes.md)

## Examples

``` r
g <- as_normalmix(krw_fit)
class(g)
#> [1] "normalmix"
# The means of the components are the grid points.
range(g$mean)
#> [1] -0.08289621  0.15852607

# The plug-in posterior distributions, with g fixed at this estimate.
plug_in <- ebnm::ebnm_npmle(krw_discrimination$yi, krw_discrimination$sei,
                            g_init = g, fix_g = TRUE)
sites <- site_summary(krw_fit)
# The posterior means are close for most firms: their absolute
# differences, in units of the posterior standard deviation.
summary(abs(coef(plug_in) - sites$mean) / sites$sd)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 6.254e-05 3.315e-03 6.467e-03 1.458e-02 1.637e-02 2.513e-01 
# The standard deviations from the fit are larger for most firms, but
# not for all.
summary(sites$sd / fitted(plug_in)$sd)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.9818  1.0034  1.0048  1.0039  1.0053  1.0724 
```
