#' Build the grid and spline basis for the effect distribution
#'
#' @description
#' The model represents the distribution of true effects, \eqn{g}, by its
#' probabilities on a grid of equally spaced points, and takes the logarithm
#' of those probabilities to be a combination of natural cubic spline
#' functions. `make_efron_grid()` builds the grid and the spline basis from
#' the estimates. [bayes_efron_fit()] calls it with its arguments
#' `grid_method`, `L`, `expansion`, `M`, `theta_true` and `bound_expansion`,
#' so the function is needed on its own only to look at a grid before
#' fitting.
#'
#' @details
#' `grid_method` chooses where the grid begins and ends. Three of the
#' methods are the rules of Lee and Sui (2025) and of the code that
#' accompanies the paper, and one is experimental.
#'
#' | `grid_method` | End points of the grid | Needs `theta_true` |
#' |:---|:---|:---|
#' | `"paper_realdata"` | The range of `theta_hat`, widened on each side by `expansion` times its width. | No |
#' | `"paper_simulation"` | The range of `theta_true`, widened by 0.5 on each side. | Yes |
#' | `"paper_sensitivity"` | The range of `theta_true`, widened on each side by `bound_expansion` times its width. | Yes |
#' | `"kl_target_experimental"` | As for `"paper_realdata"`; the number of points is set by `kappa`. | No |
#'
#' `"paper_realdata"` is the method for the analysis of data. The two
#' methods that need `theta_true`, the true effects, can be used only in a
#' simulation, where the true effects are known; they are provided so that
#' the simulations of the paper can be repeated.
#'
#' `"kl_target_experimental"` chooses the number of grid points from the
#' smallest standard error; `L` is not used. It takes the smallest number
#' for which the spacing is no more than
#' \eqn{2 \min_i(\sigma_i) \sqrt{\exp(2\kappa) - 1}}, but never fewer than
#' 51 or more than 300, so that with small standard errors and a wide range
#' of estimates the spacing is larger than the formula asks for. The
#' formula is the one that `ebnm::ebnm_scale_npmle()` uses with
#' `pointmass = FALSE`. There the distribution is a mixture of normal
#' components, each with a standard deviation of half the spacing, and
#' \eqn{\kappa}, the argument `kappa`, bounds the Kullback-Leibler
#' divergence that this approximation costs when the standard errors are
#' equal. The model of this package differs in three ways: the distribution
#' is a set of probabilities on the grid points, for which ebnm has another
#' formula; it is restricted to the family that the spline functions
#' define; and the standard errors are in general unequal. The bound has
#' not been shown to hold under these conditions, so here `kappa` only sets
#' how fine the grid is: a smaller value gives a finer grid. A message says
#' so the first time the method is used in a session. `kappa` can be set
#' only in `make_efron_grid()`; [bayes_efron_fit()] uses its default.
#'
#' The basis is `splines::ns(grid, df = M, intercept = FALSE)`. If the
#' basis and a constant column together do not have full column rank, the
#' function stops with an error; a smaller `M` avoids this.
#'
#' @param theta_hat Numeric vector of effect estimates, one for each site.
#' @param sigma Numeric vector of the standard errors of the estimates. All
#'   must be positive.
#' @param L Whole number between 51 and 300, the number of grid points, or
#'   `NULL` (the default) for 101.
#' @param expansion Number between 0 and 5. With `"paper_realdata"` and
#'   `"kl_target_experimental"`, the range of `theta_hat` is widened on each
#'   side by `expansion` times its width. Defaults to `0.5`.
#' @param kappa Number greater than 0 and less than 1, used by
#'   `"kl_target_experimental"` only: it sets the spacing that the method
#'   aims at, as described in Details. `NULL` (the default) means
#'   `1 / length(theta_hat)`.
#' @param M Whole number between 3 and 10, the number of spline functions.
#'   Defaults to `6L`.
#' @param grid_method Character string, the rule for the end points of the
#'   grid: `"paper_realdata"` (the default), `"paper_simulation"`,
#'   `"paper_sensitivity"` or `"kl_target_experimental"`. The Details of
#'   [make_efron_grid()] describe them.
#' @param theta_true Numeric vector of the true effects, of the same length
#'   as `theta_hat`. Needed by `"paper_simulation"` and
#'   `"paper_sensitivity"`, which are for simulations.
#' @param bound_expansion Number greater than 0 and at most 5. With
#'   `"paper_sensitivity"`, the range of `theta_true` is widened on each
#'   side by `bound_expansion` times its width. `NULL` (the default) means
#'   0.5.
#'
#' @return A list. `grid` holds the grid points and `B` the spline basis, a
#'   matrix with one row for each grid point and `M` columns. `L` and `M`
#'   are the numbers of grid points and spline functions, and `grid_method`
#'   is the method. `expansion` is the factor by which the range was widened
#'   (`NA` for `"paper_simulation"`), and `kappa` is `NULL` unless the method
#'   is `"kl_target_experimental"`. `attribution` is a list with the rule,
#'   `formula`, and where it comes from, `source`.
#'
#' @references
#' Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
#' deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
#' 2639. \doi{10.3390/math13162639}
#'
#' @seealso [bayes_efron_fit()], [compare_efron_grids()],
#'   `vignette("choosing-a-grid")`
#'
#' @examples
#' grid <- make_efron_grid(raudenbush1985$yi, raudenbush1985$sei)
#' range(raudenbush1985$yi)
#' range(grid$grid)
#' dim(grid$B)
#'
#' # The six spline functions on the grid.
#' matplot(grid$grid, grid$B, type = "l", lty = 1,
#'         xlab = "Effect", ylab = "Spline function")
#'
#' # A wider grid with more spline functions.
#' wide <- make_efron_grid(raudenbush1985$yi, raudenbush1985$sei,
#'                         expansion = 1, M = 8)
#' range(wide$grid)
#'
#' # In a simulation the true effects can set the range.
#' set.seed(1)
#' theta <- rnorm(50)
#' theta_hat <- rnorm(50, mean = theta, sd = 0.5)
#' sim <- make_efron_grid(theta_hat, rep(0.5, 50),
#'                        grid_method = "paper_simulation", theta_true = theta)
#' range(sim$grid)
#'
#' @export
make_efron_grid <- function(theta_hat,
                            sigma,
                            L = NULL,
                            expansion = 0.5,
                            kappa = NULL,
                            M = 6L,
                            grid_method = "paper_realdata",
                            theta_true = NULL,
                            bound_expansion = NULL) {
  grid_method <- .bef_match_grid_method(grid_method)

  theta_hat <- .bef_assert_numeric_vector(theta_hat, "theta_hat", min_len = 2L)
  sigma <- .bef_assert_numeric_vector(sigma, "sigma", min_len = 2L)
  if (length(sigma) != length(theta_hat) || any(sigma <= 0)) {
    .bef_abort_grid(
      "`sigma` must be strictly positive and the same length as `theta_hat`."
    )
  }

  L <- .bef_assert_integer_or_null(L, "L", lower = 51L, upper = 300L)
  M <- .bef_assert_integer_scalar(M, "M", lower = 3L, upper = 10L)
  expansion <- .bef_assert_number(
    expansion, "expansion", lower = 0, upper = 5
  )

  if (!is.null(kappa)) {
    kappa <- .bef_assert_number(
      kappa, "kappa", lower = 0, upper = 1, open_lower = TRUE,
      open_upper = TRUE
    )
  }
  if (!is.null(bound_expansion)) {
    bound_expansion <- .bef_assert_number(
      bound_expansion, "bound_expansion", lower = 0, upper = 5,
      open_lower = TRUE
    )
  }

  L_eff <- if (is.null(L)) 101L else L

  if (identical(grid_method, "paper_realdata")) {
    built <- .bef_grid_paper_realdata(
      theta_hat = theta_hat, L = L_eff, expansion = expansion
    )
  } else if (identical(grid_method, "paper_simulation")) {
    built <- .bef_grid_paper_simulation(
      theta_true = theta_true, theta_hat = theta_hat, L = L_eff
    )
  } else if (identical(grid_method, "paper_sensitivity")) {
    built <- .bef_grid_paper_sensitivity(
      theta_true = theta_true, theta_hat = theta_hat, L = L_eff,
      bound_expansion = bound_expansion
    )
  } else if (identical(grid_method, "kl_target_experimental")) {
    built <- .bef_grid_kl_target_experimental(
      theta_hat = theta_hat, sigma = sigma, expansion = expansion,
      kappa = kappa
    )
  } else {
    .bef_abort_grid(
      sprintf(
        "`grid_method = \"%s\"` is not implemented.",
        grid_method
      )
    )
  }

  B <- splines::ns(built$grid, df = M, intercept = FALSE)
  out <- list(
    grid = as.numeric(built$grid),
    B = B,
    M = as.integer(M),
    L = as.integer(length(built$grid)),
    expansion = built$expansion,
    kappa = built$kappa,
    grid_method = grid_method,
    attribution = built$attribution
  )

  .bef_validate_grid_return(out)
  out
}

.bef_grid_methods <- function() {
  c(
    "paper_realdata",
    "paper_simulation",
    "paper_sensitivity",
    "kl_target_experimental"
  )
}

.bef_match_grid_method <- function(grid_method) {
  if (!is.character(grid_method) || length(grid_method) != 1L ||
      is.na(grid_method)) {
    .bef_abort_grid("`grid_method` must be a single character value.")
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

.bef_grid_paper_realdata <- function(theta_hat, L, expansion) {
  range_theta <- range(theta_hat, na.rm = TRUE)
  width <- diff(range_theta)
  if (!is.finite(width) || width <= 0) {
    .bef_abort_grid("`theta_hat` must have a positive finite range.")
  }

  pad <- expansion * width
  grid <- seq(
    from = range_theta[1] - pad,
    to = range_theta[2] + pad,
    length.out = L
  )

  list(
    grid = grid,
    expansion = expansion,
    kappa = NULL,
    attribution = list(
      formula = "range(theta_hat) widened on each side by expansion * width",
      source = "Lee and Sui (2025), replication code for Section 7.2 (application) and Appendix E"
    )
  )
}

.bef_grid_paper_simulation <- function(theta_true, theta_hat, L) {
  theta_true <- .bef_require_theta_true(
    theta_true = theta_true, theta_hat = theta_hat,
    grid_method = "paper_simulation"
  )

  grid <- seq(
    from = min(theta_true) - 0.5,
    to = max(theta_true) + 0.5,
    length.out = L
  )

  list(
    grid = grid,
    expansion = NA_real_,
    kappa = NULL,
    attribution = list(
      formula = "range(theta_true) widened by 0.5 on each side",
      source = "Lee and Sui (2025), replication code for Section 5 (simulation study)"
    )
  )
}

.bef_grid_paper_sensitivity <- function(theta_true,
                                        theta_hat,
                                        L,
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
  grid <- seq(
    from = range_true[1] - range_expansion,
    to = range_true[2] + range_expansion,
    length.out = L
  )

  list(
    grid = grid,
    expansion = bound_expansion,
    kappa = NULL,
    attribution = list(
      formula = "range(theta_true) widened on each side by bound_expansion * width",
      source = "Lee and Sui (2025), Appendix C (grid sensitivity analysis)"
    )
  )
}

.bef_grid_kl_target_experimental <- function(theta_hat,
                                             sigma,
                                             expansion,
                                             kappa) {
  .bef_emit_kl_target_disclaimer()

  kappa_eff <- if (is.null(kappa)) 1 / length(theta_hat) else kappa
  range_theta <- range(theta_hat)
  width <- diff(range_theta)
  if (!is.finite(width) || width <= 0) {
    .bef_abort_grid("`theta_hat` must have a positive finite range.")
  }

  expanded_lo <- range_theta[1] - expansion * width
  expanded_hi <- range_theta[2] + expansion * width
  xrange_eff <- expanded_hi - expanded_lo
  d <- 2 * min(sigma) * sqrt(expm1(2 * kappa_eff))
  if (!is.finite(d) || d <= 0) {
    .bef_abort_grid("`kappa` and `sigma` imply a non-positive KL grid spacing.")
  }

  ideal_count <- ceiling(xrange_eff / d) + 1
  bounded_count <- as.integer(min(max(51, ideal_count), 300))
  if (!isTRUE(all.equal(ideal_count, bounded_count))) {
    .bef_inform_grid(
      sprintf(
        "kl_target_experimental grid length set to %d to enforce [51, 300] bounds.",
        bounded_count
      )
    )
  }

  grid <- seq(expanded_lo, expanded_hi, length.out = bounded_count)

  list(
    grid = grid,
    expansion = expansion,
    kappa = kappa_eff,
    attribution = list(
      formula = paste(
        "target spacing 2 * min(sigma) * sqrt(exp(2 * kappa) - 1),",
        "subject to 51 to 300 points;",
        "kappa defaults to 1 / length(theta_hat)"
      ),
      source = "ebnm::ebnm_scale_npmle() with pointmass = FALSE"
    )
  )
}

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

.bef_emit_kl_target_disclaimer <- function() {
  key <- "kl_target_experimental"
  if (isTRUE(.bayesEfron_msgs_emitted[[key]])) {
    return(invisible(FALSE))
  }

  .bef_inform_grid(.bef_kl_target_disclaimer())
  .bayesEfron_msgs_emitted[[key]] <- TRUE
  invisible(TRUE)
}

.bef_kl_target_disclaimer <- function() {
  paste(
    "`grid_method = \"kl_target_experimental\"` spaces the grid by a rule",
    "of the ebnm package, whose bound was derived for mixtures of normal",
    "components and equal standard errors. The bound has not been shown to",
    "hold for this model, and the limit of 300 grid points can make the",
    "spacing wider than the rule asks for."
  )
}

.bef_inform_grid <- function(message) {
  message(message)
  invisible(message)
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

  if (length(grid) < 51L || length(grid) > 300L) {
    .bef_abort_grid(
      sprintf(
        "Internal error in make_efron_grid(): grid length %d is outside 51 to 300 for grid_method = \"%s\".",
        length(grid), grid_method
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
