#' Construct an Efron log-spline grid
#'
#' @description
#' Build the discrete support and natural-cubic-spline basis that the
#' v0.1 random-effects Stan model uses to represent the mixing
#' distribution \eqn{g}. The function exposes the four grid recipes
#' committed by the package blueprint and the associated MDPI
#' *Mathematics* paper, together with one experimental KL-target
#' recipe.
#'
#' Most users do not need to call `make_efron_grid()` directly:
#' [bayes_efron_fit()] constructs a grid internally with the same
#' arguments. Calling this function on its own is useful when you want
#' to inspect or visualize a candidate grid before committing to a
#' fit, or to share the same grid across several fits for comparison.
#'
#' @details
#' # Grid recipes
#'
#' | Recipe | Needs `theta_true`? | Use when |
#' |:-------|:-------------------:|:---------|
#' | `"paper_realdata"` (default) | No | Real-data analysis with no oracle. Endpoints come from the observed range of `theta_hat`, expanded by `expansion`. |
#' | `"paper_simulation"` | Yes | Simulation with known truth, matched to the paper's grid rule. Endpoints come from the oracle range of `theta_true`, padded by an absolute 0.5 on each side. |
#' | `"paper_sensitivity"` | Yes | Sensitivity sweep that widens the oracle bounds by `bound_expansion`. |
#' | `"kl_target_experimental"` | No | KL-target tuning. Computes `L` from `kappa`; emits a once-per-session disclaimer because the calibration is experimental for heteroscedastic inputs. |
#'
#' Across recipes the returned grid is always strictly increasing,
#' covers the observed range, and supports a length-`M`
#' natural-cubic-spline basis (`splines::ns(grid, df = M,
#' intercept = FALSE)`).
#'
#' # Defaults and bounds
#'
#' `L` defaults to `101L` for the three paper-faithful recipes.
#' `"kl_target_experimental"` always derives the effective `L` from
#' `kappa` after validating user inputs; a supplied `L` is range-checked
#' but does not participate in the recipe's effective grid length. Both
#' `L` (51–300) and `M` (3–10) are bounded; out-of-range values are
#' rejected at the input boundary. `expansion` and `bound_expansion`
#' accept values in \eqn{[0, 5]}; `kappa` accepts values in \eqn{(0, 1)}.
#'
#' # Attribution
#'
#' Every returned grid carries an `attribution` slot recording the
#' formula and source lineage of the recipe so downstream consumers
#' can audit which paper rule produced the grid.
#'
#' @param theta_hat Numeric vector of observed effect estimates,
#'   length \eqn{\ge 2}.
#' @param sigma Numeric vector of strictly positive within-study
#'   standard errors. Same length as `theta_hat`.
#' @param L Integer grid length, or `NULL`. `NULL` uses `101L` for the
#'   three paper recipes; `"kl_target_experimental"` derives `L` from
#'   `kappa`. Bounded to \eqn{[51, 300]} when supplied.
#' @param expansion Numeric, range-relative expansion factor used by
#'   `"paper_realdata"` and `"kl_target_experimental"`. Defaults to
#'   `0.5` (50 percent expansion). Range \eqn{[0, 5]}.
#' @param kappa Numeric, KL target used by
#'   `"kl_target_experimental"`. `NULL` defaults to
#'   `1 / length(theta_hat)`. Range \eqn{(0, 1)}.
#' @param M Integer natural-cubic-spline degrees of freedom. Defaults
#'   to `6L`. Range \eqn{[3, 10]}.
#' @param grid_method Character grid recipe. One of
#'   `"paper_realdata"` (default), `"paper_simulation"`,
#'   `"paper_sensitivity"`, or `"kl_target_experimental"`.
#' @param theta_true Numeric oracle vector required by
#'   `"paper_simulation"` and `"paper_sensitivity"`; ignored by the
#'   other recipes. Same length as `theta_hat` when supplied.
#' @param bound_expansion Numeric, oracle-bound expansion factor used
#'   only by `"paper_sensitivity"`. `NULL` falls back to the recipe
#'   default of `0.5`. Range \eqn{(0, 5]}.
#'
#' @return A named list with the following fields:
#'
#'   * `grid` — numeric vector of length `L`, the discrete support of
#'     the mixing distribution.
#'   * `B` — natural-cubic-spline basis matrix from
#'     `splines::ns(grid, df = M, intercept = FALSE)`.
#'   * `M` — integer, the spline degrees of freedom actually used.
#'   * `L` — integer, the grid length actually used.
#'   * `expansion` — the range-relative expansion factor applied
#'     (recipe-dependent).
#'   * `kappa` — the KL target (only meaningful for
#'     `"kl_target_experimental"`; `NULL` otherwise).
#'   * `grid_method` — the recipe name.
#'   * `attribution` — a list recording the formula and source lineage
#'     of the recipe.
#'
#' @seealso
#'   * [bayes_efron_fit()], which constructs a grid internally and is
#'     the usual entry point.
#'   * The methodological vignette M4 (grid construction and spline
#'     basis) for the mathematical specification of each recipe.
#'
#' @examples
#' theta_hat <- c(-0.45, -0.10, 0.20, 0.55, 0.90)
#' sigma     <- c( 0.20,  0.18, 0.22, 0.16, 0.24)
#'
#' # Default real-data grid: endpoints from observed range + 50 percent expansion.
#' g_real <- make_efron_grid(theta_hat, sigma)
#' length(g_real$grid)            # 101 by default
#' dim(g_real$B)                  # L x M
#' g_real$grid_method
#'
#' # Simulation grid (oracle required): pretend the truth is known.
#' theta_true <- c(-0.40, 0.00, 0.10, 0.70, 1.00)
#' g_sim <- make_efron_grid(
#'   theta_hat   = theta_hat,
#'   sigma       = sigma,
#'   theta_true  = theta_true,
#'   grid_method = "paper_simulation"
#' )
#'
#' # Sensitivity grid: widen oracle bounds by `bound_expansion`.
#' g_sens <- make_efron_grid(
#'   theta_hat       = theta_hat,
#'   sigma           = sigma,
#'   theta_true      = theta_true,
#'   grid_method     = "paper_sensitivity",
#'   bound_expansion = 0.5
#' )
#'
#' # Experimental KL-target recipe: emits a once-per-session disclaimer.
#' g_kl <- make_efron_grid(
#'   theta_hat   = theta_hat,
#'   sigma       = sigma,
#'   grid_method = "kl_target_experimental"
#' )
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
      formula = "paper_realdata_seq",
      source = paste(
        "Lee & Sui (2025)",
        "Part_07_Real-World Application.R:228-236"
      )
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
      formula = "paper_simulation_absolute",
      source = paste(
        "Lee & Sui (2025)",
        "Part_02_Model Estimation and Performance Evaluation.R:100-106"
      )
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
      formula = "paper_sensitivity_relative",
      source = paste(
        "Lee & Sui (2025)",
        "Part_04_Grid Resolution and Bounds Sensitivity Analysis.R:112-120"
      )
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
        "d = 2 * min(sigma) * sqrt(exp(2 * kappa) - 1);",
        "kappa default = 1 / length(theta_hat)"
      ),
      source = "ebnm::ebnm_scale_npmle() (grid_selection.R:91-107)"
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
    "kl_target_experimental: ebnm KL bound assumes homoscedastic",
    "observations (`grid_selection.R:60-63`); calibration on",
    "heteroscedastic input is not validated against any published number."
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
      "make_efron_grid() postcondition violated: grid is not strictly monotone increasing."
    )
  }

  if (length(grid) < 51L || length(grid) > 300L) {
    .bef_abort_grid(
      sprintf(
        "make_efron_grid() postcondition violated: L = %d is outside [51, 300] for grid_method = \"%s\".",
        length(grid), grid_method
      )
    )
  }

  if (!is.matrix(B) || nrow(B) != length(grid) || ncol(B) != M) {
    .bef_abort_grid(
      "make_efron_grid() postcondition violated: `B` has the wrong shape."
    )
  }

  rank_aug <- qr(cbind(1, B), tol = sqrt(.Machine$double.eps))$rank
  expected_rank <- M + 1L
  if (!identical(as.integer(rank_aug), as.integer(expected_rank))) {
    rlang::abort(
      "make_efron_grid() postcondition violated: augmented spline basis is rank deficient.",
      class = c("bef_grid_rank_deficient", "bef_grid_error", "bef_error"),
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
      "make_efron_grid() postcondition violated: attribution is incomplete."
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
    class = c(class, "bef_grid_error", "bef_error")
  )
}
