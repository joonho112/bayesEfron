#' Build the grid and spline basis for the distribution of effects
#'
#' @description
#' The model represents the distribution of true effects, \eqn{g}, by its
#' probabilities on a grid of equally spaced points. The logarithm of
#' these probabilities is, up to a normalizing constant, a linear
#' combination of natural cubic spline functions (see
#' [bayes_efron_fit()]). `make_efron_grid()` builds the grid and the
#' spline basis from
#' the estimates. [bayes_efron_fit()] builds its grid in the same way from
#' its own arguments, so you only need to call this function yourself to
#' inspect a grid before fitting, or to use one grid in several fits by
#' passing it as `grid`.
#'
#' @details
#' The grid is determined by four choices.
#'
#' * The range, which is either chosen by the method named in `grid_method`
#'   or specified directly as `support`.
#' * The number of grid points, `L`.
#' * What \eqn{g} places at each grid point, `component`. With `"point"`
#'   (the default, and the model of Lee and Sui, 2025) it is a point mass.
#'   With `"normal"` it is a normal distribution whose standard deviation
#'   is `width` times the grid spacing, as in the penalized Gaussian
#'   mixture of Sui and Tipton (2026).
#' * The flexibility of \eqn{g}, which depends on the number of spline
#'   functions, `M`, and on the placement of their knots, `knots`.
#'
#' There are four methods for choosing the range. Three of them come from
#' the code that accompanies Lee and Sui (2025), and `"pilot"` narrows the
#' range of the first of them.
#'
#' | `grid_method` | Limits of the grid | Needs `theta_true` |
#' |:---|:---|:---|
#' | `"paper_realdata"` | The range of `theta_hat`, extended on each side by `expansion` times its width. | No |
#' | `"pilot"` | The range of `"paper_realdata"`, narrowed to the region where a short preliminary fit on that range (the pilot fit) places \eqn{g}. | No |
#' | `"paper_simulation"` | The range of `theta_true`, extended by 0.5 (in the units of the effects) on each side. | Yes |
#' | `"paper_sensitivity"` | The range of `theta_true`, extended on each side by `bound_expansion` times its width. | Yes |
#'
#' Use `"pilot"` (the default) or `"paper_realdata"` for real data;
#' `"paper_realdata"` does not require a model fit. The other two methods
#' need the true effects, so they can only be used in simulations. They
#' are provided so that the paper's simulations can be reproduced.
#'
#' With the default `expansion`, the range of `"paper_realdata"` is twice
#' as wide as the range of the estimates, and the estimates are themselves
#' more spread out than the true effects. With many sites, \eqn{g} then
#' occupies only a small part of the grid, most of the spline functions
#' describe the empty part, and sampling is inefficient. The `"pilot"`
#' method avoids this in two steps.
#'
#' 1. A pilot fit is run on the range of `"paper_realdata"`, with 4 chains
#'    of 500 warmup and 500 sampling iterations, 101 grid points and 6
#'    spline functions. The pilot fit always uses point masses, even when
#'    the final grid has normal components.
#' 2. The method finds the grid points that enclose the central 99.9% of
#'    the pilot estimate of \eqn{g} and extends this interval on each side
#'    by `margin` times its width, and by at least one grid spacing of the
#'    pilot fit. If an estimate lies more than 4 of its standard errors
#'    outside a limit, the limit is moved out to that estimate plus or
#'    minus 2 standard errors.
#'
#' The limits never go beyond the range of `"paper_realdata"`, so the
#' method either narrows that range or leaves it unchanged. It leaves the
#' range unchanged for [raudenbush1985] and [krw_discrimination] and
#' narrows it for [twin_towers]; see `vignette("choosing-a-grid")`. If the
#' range of `"paper_realdata"` is narrower than the median standard error,
#' no pilot fit is run. The method requires CmdStan. Specify `seed` to
#' make the result reproducible. The `pilot` component of the result
#' records what the method did.
#'
#' If `support = c(lower, upper)` is specified, the grid runs from `lower`
#' to `upper` and no method is applied, so `grid_method` must not be
#' specified. With point masses, \eqn{g} has no probability outside the
#' support, so the support is the set of values that an effect can take in
#' the model. With normal components, the support is the range of the
#' component centers, and \eqn{g} extends about 2 component standard
#' deviations further on each side. If the effects must stay within given
#' bounds, use point masses together with `reading = "step"` in the
#' summaries, because with the default `reading` a summary can extend half
#' a grid spacing beyond a limit (see [g_summary()]).
#'
#' With `knots = "grid"`, the spline basis is
#' `splines::ns(grid, df = M, intercept = FALSE)`. The interior knots are
#' at quantiles of the grid points, which makes them equally spaced over
#' the grid, and the boundary knots are at its limits. An error is thrown
#' if the basis together with a constant column does not have full column
#' rank; a smaller `M` avoids this.
#'
#' The grid must also be fine enough for the data. With point masses, the
#' likelihood of a site is evaluated only at the grid points, so a site
#' whose standard error is much smaller than the grid spacing is
#' represented poorly. In addition, a grid with only a few points across
#' the bulk of the effects cannot show the shape of \eqn{g}. The function
#' computes two ratios to check this.
#'
#' * The smallest standard error divided by the grid spacing, which should
#'   be at least 0.75.
#' * A rough estimate of the standard deviation of the effects divided by
#'   the grid spacing, which should be at least 4. With \eqn{m} the median
#'   absolute deviation of the estimates (scaled as in [stats::mad()], so
#'   that it estimates a standard deviation) and \eqn{v} the median
#'   sampling variance, this estimate is
#'   \eqn{\sqrt{m^2 - v}}{sqrt(m^2 - v)}, or \eqn{\sqrt{v}}{sqrt(v)} if
#'   \eqn{m^2 - v} is not positive.
#'
#' With `L = NULL`, the grid has 101 points if both ratios are large
#' enough, and otherwise the smallest number of points, up to 401, for
#' which they are. With normal components only the second ratio is used,
#' in an equivalent form: the standard deviation of a component should be
#' at most one sixth of the estimated standard deviation of the effects. A
#' standard error smaller than the grid spacing is not a problem in that
#' case. A warning of class `bef_warn_grid_resolution` is issued if a ratio
#' is too small for the returned grid, or if the whole grid is narrower
#' than the median standard error. The warning message suggests
#' what to change. The ratios are stored in the `resolution` component of
#' the result.
#'
#' @param theta_hat Numeric vector of effect estimates, one per site. Can
#'   also be a data frame or a `bef_data` object that contains the
#'   estimates and their standard errors (see [as_bef_data()]), in which
#'   case `sigma` is omitted.
#' @param sigma Numeric vector of standard errors of the estimates. All
#'   must be positive.
#' @param ... For `make_efron_grid()`, these dots are for future extensions
#'   and must be empty; the arguments after them must be named. The
#'   `print()` method passes them on to `format()`.
#' @param L Number of grid points, an integer from 51 to 401. If `NULL`
#'   (the default), the number is chosen from the data: 101, or more if
#'   101 points are too few for the standard errors or for the spread of
#'   the effects. See Details.
#' @param expansion Expansion factor, a number between 0 and 5. With
#'   `grid_method = "pilot"` or `"paper_realdata"`, the range of
#'   `theta_hat` is extended on each side by `expansion` times its width.
#' @param M Number of spline functions, an integer from 3 to 10. More
#'   spline functions allow \eqn{g} to have narrower features. In
#'   simulated data sets of 1,500 sites, the site intervals contained
#'   somewhat fewer of the true effects with 8 functions than with 6 (see
#'   `vignette("twin-towers")`).
#' @param grid_method Method for choosing the limits of the grid:
#'   `"pilot"` (the default), `"paper_realdata"`, `"paper_simulation"` or
#'   `"paper_sensitivity"`. See Details. `"pilot"` requires CmdStan; the
#'   others do not.
#' @param theta_true Numeric vector of true effects, the same length as
#'   `theta_hat`. Required by `"paper_simulation"` and
#'   `"paper_sensitivity"`, which are intended for simulations.
#' @param bound_expansion Expansion factor for
#'   `grid_method = "paper_sensitivity"`, a number greater than 0 and at
#'   most 5: the range of `theta_true` is extended on each side by
#'   `bound_expansion` times its width. `NULL` (the default) is equivalent
#'   to 0.5.
#' @param support Numeric vector of length 2 with the lower and upper
#'   limits of the grid. If `NULL` (the default), the limits are chosen by
#'   `grid_method`.
#' @param component Type of component at each grid point: `"point"` (the
#'   default) for a point mass, or `"normal"` for a normal distribution.
#' @param width Standard deviation of the normal components in units of the
#'   grid spacing, a number greater than 0 and at most 2. Only used with
#'   `component = "normal"`.
#' @param knots Method for placing the knots of the spline functions.
#'   Currently only `"grid"` is available; see Details.
#' @param margin Margin for `grid_method = "pilot"`, a number between 0
#'   and 1. The interval that contains the central 99.9% of the pilot
#'   estimate of \eqn{g} is extended on each side by `margin` times its
#'   width.
#' @param seed Integer seed for the pilot fit, only used with
#'   `grid_method = "pilot"`. If `NULL` (the default), a seed is taken from
#'   the clock.
#' @param x A `bef_grid` object created by `make_efron_grid()`.
#' @param digits Number of significant digits to print.
#'
#' @return An object of class `bef_grid`, which is a list with the
#'   following components.
#'   \describe{
#'     \item{`grid`}{The grid points.}
#'     \item{`B`}{The spline basis, a matrix with one row per grid point
#'       and `M` columns.}
#'     \item{`L`, `M`}{The number of grid points and the number of spline
#'       functions.}
#'     \item{`support`}{The lower and upper limits of the grid.}
#'     \item{`spacing`}{The distance between neighboring grid points.}
#'     \item{`grid_method`}{The method that chose the limits, or
#'       `"support"` if they were specified directly.}
#'     \item{`expansion`}{The factor by which the range was extended (`NA`
#'       if the method does not use one).}
#'     \item{`component`}{`"point"` or `"normal"`.}
#'     \item{`width`}{The value of `width` (`NA` for point masses).}
#'     \item{`tau_c`}{The standard deviation of a component (0 for point
#'       masses).}
#'     \item{`knots`}{A list with the method and the positions of the
#'       interior and boundary knots.}
#'     \item{`resolution`}{A list that compares the grid spacing with the
#'       data, as described in Details. It contains the two ratios,
#'       `sigma_min_over_spacing` and `scale_over_spacing`; `L_needed`,
#'       the smallest number of points for which both are large enough;
#'       `adequate`, which is `TRUE` if they are large enough for this
#'       grid; `L_chosen`, which is `TRUE` if the number of points was
#'       chosen from the data; `share_outside_estimates`, the proportion
#'       of grid points outside the range of the estimates; and
#'       `width_over_sigma`, the width of the grid divided by the median
#'       standard error.}
#'     \item{`attribution`}{A list with the formula of the method
#'       (`formula`) and its source (`source`).}
#'     \item{`pilot`}{Only for `grid_method = "pilot"`: a list that
#'       records what the method did. `wide` is the initial range and
#'       `core` is the interval that contained the central 99.9% of the
#'       pilot estimate of \eqn{g}. `narrowed` is `TRUE` if the range was
#'       narrowed, `share` is the proportion of the initial range that
#'       remains, and `protected` contains the sites whose estimates
#'       caused a limit to be moved outward. `margin` and `seed` are the
#'       values that were used. `chains`, `iter_warmup` and
#'       `iter_sampling` give the length of the pilot fit (zero if none
#'       was run), `diagnostics`
#'       contains its sampler diagnostics, and `seconds` is its run
#'       time.}
#'     \item{`data`}{A list with the estimates and standard errors from
#'       which the grid was built, used by [plot.bef_grid()]. It is not
#'       included in the grid stored in a fit.}
#'   }
#'
#' @references
#' Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
#' deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
#' 2639. \doi{10.3390/math13162639}
#'
#' Sui, D. and Tipton, E. (2026). Meta-analysis without normality:
#' estimating the true effect distribution with penalized Gaussian mixtures.
#' arXiv:2604.27887.
#'
#' @seealso [bayes_efron_fit()], [compare_efron_grids()],
#'   `vignette("choosing-a-grid")`
#'
#' @examples
#' # The default method runs the pilot fit, which needs CmdStan.
#' \dontrun{
#' make_efron_grid(twin_towers$yi, twin_towers$sei, seed = 20250627)
#' }
#'
#' # The method "paper_realdata" needs no fit.
#' grid <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata")
#' grid
#' range(krw_discrimination$yi)
#' grid$support
#' dim(grid$B)
#'
#' # The six spline functions and the knots over the estimates, in
#' # percentage points.
#' plot(grid, scale = 100)
#'
#' # A wider grid with more spline functions.
#' wide <- make_efron_grid(krw_discrimination, grid_method = "paper_realdata",
#'                         expansion = 1, M = 8)
#' wide$support
#'
#' # The limits specified directly, and a normal distribution at each grid
#' # point.
#' make_efron_grid(krw_discrimination, support = c(-0.05, 0.15),
#'                 component = "normal")
#'
#' # In a simulation the true effects can set the range.
#' sim <- make_efron_grid(twin_towers, grid_method = "paper_simulation",
#'                        theta_true = twin_towers$theta_true)
#' sim$support
#'
#' @export
make_efron_grid <- function(theta_hat,
                            sigma = NULL,
                            ...,
                            L = NULL,
                            expansion = 0.5,
                            M = 6L,
                            grid_method = "pilot",
                            theta_true = NULL,
                            bound_expansion = NULL,
                            support = NULL,
                            component = c("point", "normal"),
                            width = 2 / 3,
                            knots = "grid",
                            margin = 0.1,
                            seed = NULL) {
  .bef_grid_check_dots(list(...))
  .bef_build_grid(
    theta_hat = theta_hat, sigma = sigma, L = L, expansion = expansion,
    M = M, grid_method = grid_method, theta_true = theta_true,
    bound_expansion = bound_expansion, support = support, component = component,
    width = width, knots = knots, margin = margin, seed = seed,
    rule_given = !missing(grid_method)
  )
}

# The work of make_efron_grid(). bayes_efron_fit() calls it as well, with
# the number of chains that may run at once in the pilot fit of the rule
# "pilot". `rule_given` says whether the caller named a rule, which
# `support` does not allow. The default rule here is the one that needs no
# fit; the two exported functions pass their own.
.bef_build_grid <- function(theta_hat,
                            sigma = NULL,
                            L = NULL,
                            expansion = 0.5,
                            M = 6L,
                            grid_method = "paper_realdata",
                            theta_true = NULL,
                            bound_expansion = NULL,
                            support = NULL,
                            component = c("point", "normal"),
                            width = 2 / 3,
                            knots = "grid",
                            margin = 0.1,
                            seed = NULL,
                            rule_given = FALSE,
                            parallel_chains = NULL) {
  grid_method <- .bef_match_grid_method(grid_method)

  given <- .bef_resolve_data_args(theta_hat, sigma)
  theta_hat <- .bef_assert_numeric_vector(given$theta_hat, "theta_hat", min_len = 2L)
  sigma <- .bef_assert_numeric_vector(given$sigma, "sigma", min_len = 2L)
  if (length(sigma) != length(theta_hat) || any(sigma <= 0)) {
    .bef_abort_grid(
      "`sigma` must be strictly positive and the same length as `theta_hat`."
    )
  }

  limits <- .bef_resolution_rule()
  L <- .bef_assert_integer_or_null(L, "L", lower = limits$L_min, upper = limits$L_max)
  M <- .bef_assert_integer_scalar(M, "M", lower = 3L, upper = 10L)
  expansion <- .bef_assert_number(
    expansion, "expansion", lower = 0, upper = 5
  )

  if (!is.null(bound_expansion)) {
    bound_expansion <- .bef_assert_number(
      bound_expansion, "bound_expansion", lower = 0, upper = 5,
      open_lower = TRUE
    )
  }
  component <- .bef_match_grid_choice(component, c("point", "normal"), "component")
  knots <- .bef_match_grid_choice(knots, "grid", "knots")
  width <- .bef_assert_number(width, "width", lower = 0, upper = 2, open_lower = TRUE)
  margin <- .bef_assert_number(margin, "margin", lower = 0, upper = 1)
  seed <- .bef_assert_seed(seed)
  support <- .bef_assert_support(support)
  if (!is.null(support) && isTRUE(rule_given)) {
    .bef_abort_grid(
      "Give `support` or `grid_method`, not both: `support` sets the ends that a rule would choose."
    )
  }

  # A rule gives the two end points. The number of points follows: the one
  # given, or the one chosen from the data.
  if (!is.null(support)) {
    built <- .bef_grid_given_support(support = support)
    grid_method <- "support"
  } else if (identical(grid_method, "paper_realdata")) {
    built <- .bef_grid_paper_realdata(theta_hat = theta_hat, expansion = expansion)
  } else if (identical(grid_method, "pilot")) {
    built <- .bef_grid_pilot(
      theta_hat = theta_hat, sigma = sigma, expansion = expansion,
      margin = margin, seed = seed, parallel_chains = parallel_chains
    )
  } else if (identical(grid_method, "paper_simulation")) {
    built <- .bef_grid_paper_simulation(theta_true = theta_true, theta_hat = theta_hat)
  } else if (identical(grid_method, "paper_sensitivity")) {
    built <- .bef_grid_paper_sensitivity(
      theta_true = theta_true, theta_hat = theta_hat,
      bound_expansion = bound_expansion
    )
  } else {
    .bef_abort_grid(
      sprintf(
        "`grid_method = \"%s\"` is not implemented.",
        grid_method
      )
    )
  }

  need <- .bef_resolution_need(
    support = built$support, theta_hat = theta_hat, sigma = sigma,
    component = component, width = width
  )
  chosen <- is.null(L)
  L_eff <- if (chosen) {
    min(max(limits$L_default, need$L_needed, na.rm = TRUE), limits$L_max)
  } else {
    L
  }
  built$grid <- seq(from = built$support[1L], to = built$support[2L], length.out = L_eff)

  out <- .bef_new_grid(
    built = built, M = M, grid_method = grid_method, component = component,
    width = width, knots = knots, theta_hat = theta_hat, sigma = sigma,
    need = need, L_chosen = chosen
  )
  .bef_validate_grid_return(out)
  .bef_warn_if_coarse(out)
  out
}

# The object that make_efron_grid() returns, from the grid points that a
# rule built. Elements that version 0.3.0 had keep their names and values.
.bef_new_grid <- function(built, M, grid_method, component, width, knots,
                          theta_hat, sigma, need = NULL, L_chosen = FALSE) {
  grid <- as.numeric(built$grid)
  L <- length(grid)
  spacing <- (grid[L] - grid[1L]) / (L - 1L)
  B <- .bef_grid_basis(grid, M = M, knots = knots)
  normal <- identical(component, "normal")
  out <- c(
    list(
      grid = grid,
      B = B,
      M = as.integer(M),
      L = as.integer(L),
      expansion = built$expansion,
      grid_method = grid_method,
      attribution = built$attribution,
      support = c(grid[1L], grid[L]),
      spacing = spacing,
      component = component,
      width = if (normal) width else NA_real_,
      tau_c = if (normal) width * spacing else 0,
      knots = list(
        rule = knots,
        interior = as.numeric(attr(B, "knots")),
        boundary = as.numeric(attr(B, "Boundary.knots"))
      ),
      resolution = .bef_grid_resolution(
        grid, spacing, theta_hat, sigma, component = component, width = width,
        need = need, L_chosen = L_chosen
      )
    ),
    # What the rule "pilot" did; other rules have no such element.
    if (!is.null(built$pilot)) list(pilot = built$pilot),
    # The estimates that the grid was built from, for its plot.
    list(data = list(theta_hat = as.numeric(theta_hat), sigma = as.numeric(sigma)))
  )
  structure(out, class = "bef_grid")
}

# The spline basis on the grid. Rule "grid" is that of version 0.3.0: the
# matrix depends on the number of grid points and of functions only, and not
# on where the grid lies.
.bef_grid_basis <- function(grid, M, knots = "grid") {
  splines::ns(grid, df = M, intercept = FALSE)
}

# The numbers of the rule for the spacing. With point masses the smallest
# standard error is at least `sigma` spacings and the scale of the effects
# at least `scale` spacings; with normal components of the default width the
# second condition alone. A grid narrower than `width` median standard
# errors is reported. The number of points is `L_default` unless more are
# needed, and never more than `L_max`.
.bef_resolution_rule <- function() {
  list(sigma = 0.75, scale = 4, width = 1, L_min = 51L, L_default = 101L, L_max = 401L)
}

# The largest spacing that the rule allows for these data, and the smallest
# number of points that gives it on `support`.
.bef_resolution_need <- function(support, theta_hat, sigma, component, width) {
  rule <- .bef_resolution_rule()
  scale <- .bef_robust_scale(theta_hat, sigma)
  spacing_max <- if (identical(component, "normal")) {
    # A component of standard deviation `width` spacings is at most a sixth
    # of the scale: at the default width, 2 / 3, the scale is four spacings.
    scale * (2 / 3) / (rule$scale * width)
  } else {
    min(min(sigma) / rule$sigma, scale / rule$scale)
  }
  points <- (support[2L] - support[1L]) / spacing_max
  list(
    spacing_max = spacing_max,
    # The tolerance keeps a ratio that is a whole number from rounding up.
    L_needed = if (is.finite(points) && points < .Machine$integer.max / 2) {
      as.integer(ceiling(points * (1 - 1e-12))) + 1L
    } else {
      NA_integer_
    }
  )
}

# How the spacing of a grid compares with the data.
.bef_grid_resolution <- function(grid, spacing, theta_hat, sigma, component = "point",
                                 width = NA_real_, need = NULL, L_chosen = FALSE) {
  rule <- .bef_resolution_rule()
  scale <- .bef_robust_scale(theta_hat, sigma)
  estimates <- range(theta_hat)
  L <- length(grid)
  if (is.null(need)) {
    need <- .bef_resolution_need(c(grid[1L], grid[L]), theta_hat, sigma, component, width)
  }
  list(
    sigma_min_over_spacing = min(sigma) / spacing,
    scale_over_spacing = scale / spacing,
    share_outside_estimates = mean(grid < estimates[1L] | grid > estimates[2L]),
    width_over_sigma = (grid[L] - grid[1L]) / stats::median(sigma),
    L_chosen = isTRUE(L_chosen),
    L_needed = need$L_needed,
    # The spacing itself is compared, with a tolerance for rounding.
    adequate = spacing <= need$spacing_max * (1 + 1e-9)
  )
}

# One warning for a grid that is too coarse or too narrow for the data. It
# names every ratio that is too small and says what would help.
.bef_warn_if_coarse <- function(grid) {
  rule <- .bef_resolution_rule()
  res <- grid$resolution
  normal <- identical(grid$component, "normal")
  number <- function(x) format(signif(x, 2), trim = TRUE)
  problems <- character()
  if (!isTRUE(res$adequate)) {
    if (!normal && res$sigma_min_over_spacing < rule$sigma) {
      problems <- c(problems, sprintf(
        "The smallest standard error is %s spacings of the grid; below %s the likelihood of such a site falls between the grid points.",
        number(res$sigma_min_over_spacing), rule$sigma
      ))
    }
    if (!normal && res$scale_over_spacing < rule$scale) {
      problems <- c(problems, sprintf(
        "The spread of the effects is about %s spacings of the grid; below %s the grid is too coarse for the shape of their distribution.",
        number(res$scale_over_spacing), rule$scale
      ))
    }
    if (normal) {
      problems <- c(problems, sprintf(
        "A normal component is %s times as wide as the spread of the effects; above %s the components blur their distribution.",
        number(grid$width / res$scale_over_spacing), number((2 / 3) / rule$scale)
      ))
    }
    problems <- c(problems, if (is.na(res$L_needed) || res$L_needed > rule$L_max) {
      if (normal) {
        sprintf("No `L` up to %d is enough; a smaller `width` or a narrower range (`support`) would help.",
                rule$L_max)
      } else {
        sprintf("No `L` up to %d is enough. `component = \"normal\"` puts a normal distribution at each grid point, for which a small standard error does no harm, and `support` can set a narrower range.",
                rule$L_max)
      }
    } else {
      sprintf("`L = %d` is enough, and `L = NULL` chooses it.", res$L_needed)
    })
  }
  if (res$width_over_sigma < rule$width) {
    problems <- c(problems, sprintf(
      "The grid is %s wide, less than the median standard error (%s): the estimates agree more closely than their standard errors can explain, and the data say little about a distribution on so narrow a range. Check the estimates and the standard errors, or give a range in `support`.",
      number(grid$support[2L] - grid$support[1L]),
      number((grid$support[2L] - grid$support[1L]) / res$width_over_sigma)
    ))
  }
  if (length(problems) == 0L) {
    return(invisible(FALSE))
  }
  .bef_warn(
    paste(problems, collapse = " "),
    "bef_warn_grid_resolution",
    L = grid$L,
    L_needed = res$L_needed,
    sigma_min_over_spacing = res$sigma_min_over_spacing,
    scale_over_spacing = res$scale_over_spacing,
    width_over_sigma = res$width_over_sigma
  )
  invisible(TRUE)
}

# A robust estimate of the standard deviation of the true effects: the
# square root of the squared median absolute deviation of the estimates less
# the median sampling variance, or, where that is not positive, of the median
# sampling variance (Sui and Tipton 2026, Section 2.2).
.bef_robust_scale <- function(theta_hat, sigma) {
  typical <- stats::median(sigma^2)
  excess <- stats::mad(theta_hat)^2 - typical
  sqrt(if (is.finite(excess) && excess > 0) excess else typical)
}

.bef_grid_given_support <- function(support) {
  list(
    support = support,
    expansion = NA_real_,
    attribution = list(
      formula = "ends given in `support`",
      source = "the call"
    )
  )
}

.bef_assert_support <- function(support) {
  if (is.null(support)) {
    return(NULL)
  }
  if (!is.numeric(support) || length(support) != 2L || any(!is.finite(support)) ||
      support[1L] >= support[2L]) {
    .bef_abort_grid(
      "`support` must be two finite numbers, the lower end of the grid and then the upper end."
    )
  }
  as.numeric(support)
}

.bef_match_grid_choice <- function(value, choices, arg) {
  if (identical(value, choices)) {
    return(choices[[1L]])
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    .bef_abort_grid(
      sprintf("`%s` must be %s.", arg,
              paste(sprintf('"%s"', choices), collapse = " or "))
    )
  }
  value
}

# Arguments after the data must be named in full. `grid` has a message of
# its own: a grid of the user's own points is what a caller is most likely
# to ask for that the function does not take.
.bef_grid_check_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }
  given <- names(dots)
  if (!is.null(given) && "kappa" %in% given) {
    .bef_abort_grid(.bef_kl_target_removed())
  }
  if (!is.null(given) && "grid" %in% given) {
    .bef_abort_grid(
      "`grid` is not an argument: the grid points are equally spaced. Give the two ends in `support` and the number of points in `L`."
    )
  }
  named <- if (is.null(given)) character() else given[nzchar(given)]
  .bef_abort_grid(
    if (length(named) > 0L) {
      sprintf("Unused arguments: %s.", paste(sprintf("`%s`", named), collapse = ", "))
    } else {
      "The settings of the grid must be given by name, as in `L = 101`."
    }
  )
}

#' @describeIn make_efron_grid Returns the lines that `print()` shows as a
#'   character vector.
#' @export
format.bef_grid <- function(x, digits = 3L, ...) {
  # Each number on its own, so that one does not set the decimals of another.
  number <- function(value) {
    vapply(signif(value, digits), format, character(1L), trim = TRUE)
  }
  rule <- if (identical(x$grid_method, "support")) {
    "ends given in `support`"
  } else if (identical(x$grid_method, "pilot")) {
    sprintf(
      "\"pilot\", expansion %s, margin %s: %s", number(x$expansion), number(x$pilot$margin),
      if (isTRUE(x$pilot$narrowed)) {
        sprintf("narrowed to %s%% of the range of \"paper_realdata\"",
                number(100 * x$pilot$share))
      } else {
        "not narrowed, the range of \"paper_realdata\""
      }
    )
  } else {
    sprintf("\"%s\"%s", x$grid_method,
            if (isTRUE(is.finite(x$expansion))) sprintf(", expansion %s", number(x$expansion)) else "")
  }
  component <- if (identical(x$component, "normal")) {
    sprintf("normal, standard deviation %s (%s of the spacing)", number(x$tau_c), number(x$width))
  } else {
    "point masses"
  }
  c(
    sprintf("<bef_grid> %d points from %s to %s, spacing %s",
            x$L, number(x$support[1L]), number(x$support[2L]), number(x$spacing)),
    sprintf("  Rule:       %s", rule),
    sprintf("  Component:  %s", component),
    sprintf("  Basis:      %d natural cubic spline functions, knots by rule \"%s\"",
            x$M, x$knots$rule),
    sprintf("  Interior knots: %s", paste(number(x$knots$interior), collapse = ", ")),
    sprintf("  Resolution: smallest standard error %s spacings, scale of the effects %s spacings%s",
            number(x$resolution$sigma_min_over_spacing), number(x$resolution$scale_over_spacing),
            if (isFALSE(x$resolution$adequate)) " (too coarse)" else ""),
    sprintf("  %s%% of the grid points lie outside the range of the estimates",
            number(100 * x$resolution$share_outside_estimates))
  )
}

#' @describeIn make_efron_grid Prints the limits and the spacing of the
#'   grid, the method that chose the limits, the components, the spline
#'   functions with their knots, and how the spacing compares with the
#'   data.
#' @export
print.bef_grid <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

.bef_grid_methods <- function() {
  c(
    "paper_realdata",
    "pilot",
    "paper_simulation",
    "paper_sensitivity"
  )
}

# What is said about the rule that version 0.4 removed.
.bef_kl_target_removed <- function() {
  paste(
    "`grid_method = \"kl_target_experimental\"` and its argument `kappa` were removed in",
    "version 0.4: the formula was made for normal components and equal standard errors.",
    "Leave `L` out, and the number of grid points is chosen from the smallest standard",
    "error and the spread of the effects."
  )
}

.bef_match_grid_method <- function(grid_method) {
  if (!is.character(grid_method) || length(grid_method) != 1L ||
      is.na(grid_method)) {
    .bef_abort_grid("`grid_method` must be a single character value.")
  }
  if (identical(grid_method, "kl_target_experimental")) {
    .bef_abort_grid(.bef_kl_target_removed())
  }
  if (!grid_method %in% .bef_grid_methods()) {
    .bef_abort_grid(
      sprintf(
        "`grid_method` must be one of: %s.",
        paste(sprintf('"%s"', .bef_grid_methods()), collapse = ", ")
      )
    )
  }
  grid_method
}

.bef_grid_paper_realdata <- function(theta_hat, expansion) {
  list(
    support = .bef_realdata_range(theta_hat, expansion),
    expansion = expansion,
    attribution = list(
      formula = "range(theta_hat) widened on each side by expansion * width",
      source = "Lee and Sui (2025), replication code for Section 7.2 (application) and Appendix E"
    )
  )
}

.bef_grid_paper_simulation <- function(theta_true, theta_hat) {
  theta_true <- .bef_require_theta_true(
    theta_true = theta_true, theta_hat = theta_hat,
    grid_method = "paper_simulation"
  )

  list(
    support = c(min(theta_true) - 0.5, max(theta_true) + 0.5),
    expansion = NA_real_,
    attribution = list(
      formula = "range(theta_true) widened by 0.5 on each side",
      source = "Lee and Sui (2025), replication code for Section 5 (simulation study)"
    )
  )
}

.bef_grid_paper_sensitivity <- function(theta_true,
                                        theta_hat,
                                        bound_expansion) {
  theta_true <- .bef_require_theta_true(
    theta_true = theta_true, theta_hat = theta_hat,
    grid_method = "paper_sensitivity"
  )
  bound_expansion <- if (is.null(bound_expansion)) 0.5 else bound_expansion

  range_true <- range(theta_true)
  width <- diff(range_true)
  if (!is.finite(width) || width <= 0) {
    .bef_abort_grid("`theta_true` must have a positive finite range.")
  }
  range_expansion <- width * bound_expansion

  list(
    support = c(range_true[1] - range_expansion, range_true[2] + range_expansion),
    expansion = bound_expansion,
    attribution = list(
      formula = "range(theta_true) widened on each side by bound_expansion * width",
      source = "Lee and Sui (2025), Appendix C (grid sensitivity analysis)"
    )
  )
}

# What was said once in the session.
.bayesEfron_msgs_emitted <- new.env(parent = emptyenv())

.bef_require_theta_true <- function(theta_true, theta_hat, grid_method) {
  if (is.null(theta_true)) {
    .bef_abort_grid(
      sprintf(
        "`theta_true` is required for `grid_method = \"%s\"`.",
        grid_method
      ),
      class = "bef_err_grid_oracle_required"
    )
  }

  theta_true <- .bef_assert_numeric_vector(
    theta_true, "theta_true", min_len = 1L
  )
  if (length(theta_true) != length(theta_hat)) {
    .bef_abort_grid("`theta_true` must be the same length as `theta_hat`.")
  }
  theta_true
}

.bef_validate_grid_return <- function(out) {
  grid <- out$grid
  B <- out$B
  M <- out$M
  grid_method <- out$grid_method
  attribution <- out$attribution

  if (!is.numeric(grid) || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    .bef_abort_grid(
      "Internal error in make_efron_grid(): the grid is not strictly increasing."
    )
  }

  limits <- .bef_resolution_rule()
  if (length(grid) < limits$L_min || length(grid) > limits$L_max) {
    .bef_abort_grid(
      sprintf(
        "Internal error in make_efron_grid(): grid length %d is outside %d to %d for grid_method = \"%s\".",
        length(grid), limits$L_min, limits$L_max, grid_method
      )
    )
  }

  if (!is.matrix(B) || nrow(B) != length(grid) || ncol(B) != M) {
    .bef_abort_grid(
      "Internal error in make_efron_grid(): the basis matrix has the wrong dimensions."
    )
  }

  rank_aug <- qr(cbind(1, B), tol = sqrt(.Machine$double.eps))$rank
  expected_rank <- M + 1L
  if (!identical(as.integer(rank_aug), as.integer(expected_rank))) {
    rlang::abort(
      "The spline basis does not have full column rank on this grid; use a smaller `M`.",
      class = c("bef_grid_rank_deficient", "bef_grid_error", "bef_error"),
      call = .bef_user_frame(),
      rank = as.integer(rank_aug),
      expected_rank = as.integer(expected_rank),
      L = length(grid),
      grid_range = range(grid),
      grid_method = grid_method
    )
  }

  if (!is.list(attribution) ||
      !is.character(attribution$formula) ||
      length(attribution$formula) != 1L ||
      !nzchar(attribution$formula) ||
      !is.character(attribution$source) ||
      length(attribution$source) != 1L ||
      !nzchar(attribution$source)) {
    .bef_abort_grid(
      "Internal error in make_efron_grid(): the attribution entry is incomplete."
    )
  }

  invisible(out)
}

.bef_assert_numeric_vector <- function(x, arg, min_len = 1L) {
  if (!is.numeric(x) || length(x) < min_len || any(!is.finite(x))) {
    .bef_abort_grid(
      sprintf("`%s` must be a finite numeric vector of length at least %d.", arg, min_len)
    )
  }
  as.numeric(x)
}

# A seed: `NULL`, or a whole number that CmdStan accepts.
.bef_assert_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    .bef_abort_grid("`seed` must be a single whole number that is not negative, or `NULL`.")
  }
  as.integer(seed)
}

.bef_assert_integer_or_null <- function(x, arg, lower, upper) {
  if (is.null(x)) {
    return(NULL)
  }
  .bef_assert_integer_scalar(x, arg, lower = lower, upper = upper)
}

.bef_assert_integer_scalar <- function(x, arg, lower, upper) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x)) {
    .bef_abort_grid(sprintf("`%s` must be a single integer.", arg))
  }

  x <- as.integer(x)
  if (x < lower || x > upper) {
    .bef_abort_grid(
      sprintf("`%s` must be between %d and %d.", arg, lower, upper)
    )
  }
  x
}

.bef_assert_number <- function(x,
                               arg,
                               lower,
                               upper,
                               open_lower = FALSE,
                               open_upper = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    .bef_abort_grid(sprintf("`%s` must be a single finite number.", arg))
  }

  lower_ok <- if (open_lower) x > lower else x >= lower
  upper_ok <- if (open_upper) x < upper else x <= upper
  if (!lower_ok || !upper_ok) {
    left <- if (open_lower) "(" else "["
    right <- if (open_upper) ")" else "]"
    .bef_abort_grid(
      sprintf("`%s` must be in %s%s, %s%s.", arg, left, lower, upper, right)
    )
  }
  as.numeric(x)
}

.bef_abort_grid <- function(message, class = "bef_invalid_args") {
  rlang::abort(
    message,
    class = c(class, "bef_grid_error", "bef_error"),
    call = .bef_user_frame()
  )
}
