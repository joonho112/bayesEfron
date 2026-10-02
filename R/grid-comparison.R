# compare_efron_grids(): fit several grid specifications and rank them.

#' Compare grids and numbers of spline functions by cross-validation
#'
#' @description
#' The number of grid points, the width of the grid and the number of
#' spline functions are choices of the analyst.
#' `compare_efron_grids()` fits the model once for every combination of the
#' values given in `L`, `M` and `expansion` and ranks the fits by
#' `elpd_loo`, the expected log predictive density for a new site as
#' estimated by leave-one-out cross-validation ([loo.bef_fit()]).
#'
#' @details
#' The table is ordered from the largest `elpd_loo` to the smallest.
#' `elpd_diff` is the difference from the first row and `se_diff` is its
#' standard error. A difference of less than about two standard errors
#' does not separate two specifications. With few sites the standard error
#' is itself uncertain, and the rule is a rough guide.
#'
#' A fit that has not converged should not be chosen for its `elpd_loo`.
#' Each row therefore carries the largest R-hat, the smallest bulk
#' effective sample size, the number of divergent transitions and the
#' largest Pareto k of its fit, and `print()` names the rows with an R-hat
#' above 1.05, a divergent transition or a Pareto k above 0.7.
#'
#' With the default arguments nine models are fitted, each with the
#' sampler settings of [bayes_efron_fit()] unless they are changed through
#' `...`. A fit that fails gives a warning and a row of missing values.
#' The function needs the loo package.
#'
#' @inheritParams bayes_efron_fit
#' @param L Vector of whole numbers between 51 and 300, the numbers of grid
#'   points to try. Defaults to `c(51, 101, 201)`.
#' @param M Vector of whole numbers between 3 and 10, the numbers of spline
#'   functions to try. Defaults to `c(4, 6, 8)`.
#' @param expansion Numeric vector of values between 0 and 5, the
#'   expansion factors to try. Defaults to `0.5`.
#' @param ... For `compare_efron_grids()`, other arguments of
#'   [bayes_efron_fit()], such as `chains`, `iter_sampling` and `seed`. Not
#'   used by `print()`.
#' @param x A `bef_grid_comparison` object.
#'
#' @return A data frame of class `bef_grid_comparison` with one row for
#'   each specification, the best first. Its columns are `L`, `M` and
#'   `expansion`; `elpd_loo`, its standard error `se_elpd_loo` and `p_loo`,
#'   the effective number of parameters as the loo package defines it,
#'   which is not the `effective_params` of a fit; `elpd_diff` and `se_diff`;
#'   `max_rhat`, `min_ess_bulk`, `divergences` and `max_pareto_k`; and
#'   `sampler_seconds`, the time that sampling took.
#'
#' @seealso [loo.bef_fit()], [make_efron_grid()],
#'   `vignette("choosing-a-grid")`
#'
#' @examples
#' # Nine fits of the teacher expectancy data. It needs CmdStan.
#' \dontrun{
#' compare_efron_grids(raudenbush1985$yi, raudenbush1985$sei, seed = 1985)
#' }
#'
#' @export
compare_efron_grids <- function(theta_hat,
                                sigma,
                                L = c(51L, 101L, 201L),
                                M = c(4L, 6L, 8L),
                                expansion = 0.5,
                                grid_method = "paper_realdata",
                                ...) {
  .bef_require_loo()
  L <- vapply(L, function(x) .bef_assert_integer_scalar(x, "L", lower = 51L, upper = 300L),
              integer(1))
  M <- vapply(M, function(x) .bef_assert_integer_scalar(x, "M", lower = 3L, upper = 10L),
              integer(1))
  expansion <- vapply(expansion, function(x) .bef_assert_number(x, "expansion", lower = 0, upper = 5),
                      numeric(1))

  specs <- expand.grid(L = L, M = M, expansion = expansion,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  rows <- vector("list", nrow(specs))
  loos <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    fit <- tryCatch(
      suppressWarnings(bayes_efron_fit(
        theta_hat = theta_hat, sigma = sigma,
        L = spec$L, M = spec$M, expansion = spec$expansion,
        grid_method = grid_method, ...
      )),
      error = function(err) {
        .bef_warn(
          sprintf("Fit failed for L = %d, M = %d: %s",
                  spec$L, spec$M, conditionMessage(err)),
          "bef_grid_comparison_fit_failed"
        )
        NULL
      }
    )
    if (is.null(fit)) {
      rows[[i]] <- data.frame(
        L = spec$L, M = spec$M, expansion = spec$expansion,
        elpd_loo = NA_real_, se_elpd_loo = NA_real_, p_loo = NA_real_,
        max_rhat = NA_real_, min_ess_bulk = NA_real_,
        divergences = NA_real_, max_pareto_k = NA_real_,
        sampler_seconds = NA_real_
      )
      next
    }

    lo <- loo::loo(fit)
    d <- diagnose(fit)
    loos[[i]] <- lo
    rows[[i]] <- data.frame(
      L = spec$L, M = spec$M, expansion = spec$expansion,
      elpd_loo = lo$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = lo$estimates["elpd_loo", "SE"],
      p_loo = lo$estimates["p_loo", "Estimate"],
      max_rhat = suppressWarnings(max(d$rhat, na.rm = TRUE)),
      min_ess_bulk = suppressWarnings(min(d$ess_bulk, na.rm = TRUE)),
      divergences = sum(d$divergences, na.rm = TRUE),
      max_pareto_k = suppressWarnings(max(lo$diagnostics$pareto_k, na.rm = TRUE)),
      sampler_seconds = fit$metadata$runtime_seconds
    )
    rm(fit)
    invisible(gc(verbose = FALSE))
  }

  out <- do.call(rbind, rows)
  ok <- which(!is.na(out$elpd_loo))
  out$elpd_diff <- NA_real_
  out$se_diff <- NA_real_
  if (length(ok) > 0L) {
    differences <- .bef_elpd_differences(loos[ok])
    out$elpd_diff[ok] <- differences$elpd_diff
    out$se_diff[ok] <- differences$se_diff
  }

  out <- out[order(-out$elpd_loo, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out <- out[, c("L", "M", "expansion", "elpd_loo", "se_elpd_loo", "p_loo",
                 "elpd_diff", "se_diff", "max_rhat", "min_ess_bulk",
                 "divergences", "max_pareto_k", "sampler_seconds")]
  structure(out, class = c("bef_grid_comparison", "data.frame"))
}

# The difference in elpd_loo between each fit and the best one, with its
# standard error, from the pointwise values: the sum of the differences and
# sqrt(n) times their standard deviation, as loo::loo_compare() computes them.
# The values are returned in the order of `loos`, so that they do not depend
# on how loo_compare() labels its rows, which has changed between versions of
# loo.
.bef_elpd_differences <- function(loos) {
  pointwise <- lapply(loos, function(x) x$pointwise[, "elpd_loo"])
  totals <- vapply(pointwise, sum, numeric(1))
  reference <- pointwise[[which.max(totals)]]
  differences <- lapply(pointwise, function(x) x - reference)
  list(
    elpd_diff = vapply(differences, sum, numeric(1)),
    se_diff = vapply(
      differences,
      function(d) sqrt(length(d)) * stats::sd(d),
      numeric(1)
    )
  )
}

#' @rdname compare_efron_grids
#' @export
print.bef_grid_comparison <- function(x, ...) {
  cat("<bef_grid_comparison>", nrow(x), "specifications, best first\n\n")
  show <- as.data.frame(x)
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], function(v) round(v, 3))
  print(show, row.names = FALSE)

  flagged <- which(
    (!is.na(x$max_rhat) & x$max_rhat > 1.05) |
      (!is.na(x$divergences) & x$divergences > 0) |
      (!is.na(x$max_pareto_k) & x$max_pareto_k > 0.7)
  )
  if (length(flagged) > 0L) {
    cat("\nRows with convergence or PSIS problems:",
        paste(flagged, collapse = ", "),
        "\nA specification that did not converge should not be selected on elpd_loo alone.\n")
  }
  invisible(x)
}
