# compare_efron_grids(): fit several grid specifications and set them side by
# side.

#' Compare fits on several grids (sensitivity analysis)
#'
#' @description
#' `compare_efron_grids()` fits the model on several grids and returns a
#' table that compares the fits. The range of the grid, the number of grid
#' points, the type of component and the number of spline functions are
#' choices of the analyst, and the table shows how sensitive the results
#' are to these choices. The model is fitted once for each combination of
#' the specified values, called a specification. For each specification,
#' the table contains the posterior mean of the standard deviation of
#' \eqn{g}, the posterior means of a few sites, sampler and grid
#' diagnostics, and `elpd_loo`, the expected log predictive density for a
#' new site estimated by leave-one-out cross-validation (see
#' [loo.bef_fit()]).
#'
#' @details
#' The column `sd_g` and the site columns show how the estimate of
#' \eqn{g} and the site effects change with the grid. By default the site
#' columns are for the sites with the smallest and the largest estimates,
#' which are the sites that a change of the grid affects most. To vary
#' only how far the grid extends beyond the estimates, specify several
#' values of `expansion` and a single value each for `L` and `M`.
#'
#' The rows are sorted by `elpd_loo`, with the largest value first.
#' `elpd_diff` is the difference in `elpd_loo` from the first row, and
#' `se_diff` is the standard error of this difference. A difference of
#' less than about 2 standard errors is not enough to distinguish two
#' specifications, which is the usual case when there are few sites. Note
#' that `elpd_loo` measures how well a fit predicts the estimate of a new
#' site, and this is hardly affected by the limits of the grid. For
#' example, in the comparison of three ranges for [raudenbush1985] that is
#' loaded in the examples, no fit differs from the best by 2 standard
#' errors, although the posterior mean of the study with the largest
#' estimate (column `site_4`) is much lower with the widest range than
#' with the narrowest. `print()` adds a note if no specification can be
#' distinguished from the best one.
#'
#' A fit that has not converged should not be selected for its `elpd_loo`.
#' Each row therefore includes the largest R-hat, the smallest bulk
#' effective sample size, the number of divergent transitions and the
#' largest Pareto k of its fit. `print()` lists the rows with an R-hat
#' above 1.05, a divergent transition or a Pareto k above 0.7.
#'
#' Two more columns check the grid of each fit. `end_mass` is the largest
#' posterior mass that any site has in the outer 2% of the range at one
#' end of the grid, that is, at its lower or upper limit, as reported by
#' [diagnose()]. A value above 0.05 means that the result for that site
#' depends on the limits of the grid. `adequate` is `TRUE` if the grid is
#' fine enough for the data, as defined in the Details of
#' [make_efron_grid()]. Warnings from the individual fits are suppressed.
#'
#' With the default arguments, 9 models are fitted. Each fit uses the
#' default sampler settings of [bayes_efron_fit()] unless others are
#' specified through `...`. A message is shown as each fit starts (in
#' these messages, "rule" means the grid method). If a fit fails, a
#' warning is issued and its row has missing values, but the remaining
#' fits are still run. The function requires the loo package.
#'
#' @inheritParams bayes_efron_fit
#' @param L Numbers of grid points to compare, a vector of integers from 51
#'   to 401.
#' @param M Numbers of spline functions to compare, a vector of integers
#'   from 3 to 10.
#' @param expansion Expansion factors to compare, a numeric vector with
#'   values between 0 and 5. See the `expansion` argument of
#'   [make_efron_grid()].
#' @param grid_method Methods for choosing the limits of the grid, a
#'   character vector with one or more of the methods listed in the
#'   Details of [make_efron_grid()]. The default is `"paper_realdata"`, not
#'   `"pilot"` as in [bayes_efron_fit()], so that `expansion` determines
#'   the range directly.
#' @param ... Further arguments passed to [bayes_efron_fit()] for every
#'   fit, e.g., `chains`, `iter_sampling` or `seed`. Not used by the
#'   `print()` method.
#' @param support Numeric vector of length 2 with the lower and upper
#'   limits of the grid, or a list of such vectors to compare several
#'   ranges. If `support` is specified, `grid_method` and `expansion` are
#'   not used.
#' @param component One or both of `"point"` (the default) and `"normal"`,
#'   the types of component to compare.
#' @param sites Sites whose posterior means are reported, specified as a
#'   numeric vector of positions in the data or as a character vector of
#'   site labels. If `NULL` (the default), the sites with the smallest and
#'   the largest estimates are used.
#' @param keep_fits If `TRUE`, the fitted objects are stored in the `fits`
#'   attribute of the result.
#' @param x A `bef_grid_comparison` object returned by
#'   `compare_efron_grids()`.
#'
#' @return A data frame of class `bef_grid_comparison` with one row per
#'   specification, sorted by decreasing `elpd_loo` (best first), and the
#'   following columns.
#'   \describe{
#'     \item{`spec`}{The number of the specification, in the order in which
#'       the specifications were fitted.}
#'     \item{`L`, `M`, `expansion`, `grid_method`, `component`}{The
#'       settings of the specification.}
#'     \item{`lower`, `upper`}{The lower and upper limits of the grid.}
#'     \item{`elpd_loo`}{The expected log predictive density, estimated by
#'       leave-one-out cross-validation.}
#'     \item{`se_elpd_loo`}{The standard error of `elpd_loo`.}
#'     \item{`p_loo`}{The effective number of parameters as defined by the
#'       loo package. This is not the `effective_params` of a fit.}
#'     \item{`elpd_diff`}{The difference in `elpd_loo` from the first row.}
#'     \item{`se_diff`}{The standard error of `elpd_diff`.}
#'     \item{`mean_g`}{The posterior mean of the mean of \eqn{g}.}
#'     \item{`sd_g`}{The posterior mean of the standard deviation of
#'       \eqn{g}.}
#'     \item{\code{site_<label>}}{The posterior mean of the effect of a
#'       selected site, one column per site, named after the site label (e.g.,
#'       `site_4`).}
#'     \item{`max_rhat`}{The largest R-hat of the fit.}
#'     \item{`min_ess_bulk`}{The smallest bulk effective sample size.}
#'     \item{`divergences`}{The number of divergent transitions.}
#'     \item{`max_pareto_k`}{The largest Pareto k.}
#'     \item{`end_mass`}{The largest posterior mass that any site has in
#'       the outer 2% of the range at one end of the grid (see Details).}
#'     \item{`adequate`}{`TRUE` if the grid is fine enough for the data
#'       (see Details).}
#'     \item{`sampler_seconds`}{The sampling time in seconds.}
#'   }
#'   The attribute `g` is a list, named by `spec`, with the density of
#'   \eqn{g} for each specification as returned by [g_summary()]. The
#'   attribute `sites` is a data frame with the rows of [site_summary()]
#'   for the selected sites and a `spec` column. With `keep_fits = TRUE`,
#'   the attribute `fits` contains the fitted objects.
#'
#' @seealso [plot.bef_grid_comparison()], [loo.bef_fit()],
#'   [make_efron_grid()], [diagnose()], `vignette("choosing-a-grid")`
#'
#' @examples
#' # Both need CmdStan.
#' \dontrun{
#' # Nine fits of the firm-level callback gaps.
#' compare_efron_grids(krw_discrimination, seed = 20250807)
#'
#' # How far the grid extends beyond the estimates, for the 19 studies.
#' compare_efron_grids(raudenbush1985, L = 101, M = 6,
#'                     expansion = c(0.25, 0.5, 1), sites = c(4, 8),
#'                     seed = 1985)
#' }
#'
#' # The result of the second call is stored with the package.
#' comparison <- readRDS(system.file("examples", "raudenbush_ranges.rds",
#'                                   package = "bayesEfron"))
#' comparison
#' plot(comparison)
#'
#' @export
compare_efron_grids <- function(theta_hat,
                                sigma = NULL,
                                L = c(51L, 101L, 201L),
                                M = c(4L, 6L, 8L),
                                expansion = 0.5,
                                grid_method = "paper_realdata",
                                ...,
                                support = NULL,
                                component = "point",
                                sites = NULL,
                                keep_fits = FALSE) {
  .bef_require_loo()
  # `grid` would otherwise be taken for `grid_method`, whose name it begins.
  if ("grid" %in% names(sys.call())) {
    .bef_abort_invalid_args(
      "`grid` cannot be given here: the function builds the grids that it compares.",
      arg = "grid"
    )
  }
  # What every fit would reject is rejected once, before any fit is made.
  given <- .bef_resolve_data_args(theta_hat, sigma)
  theta_hat <- .bef_validate_theta_hat(given$theta_hat)
  sigma <- .bef_validate_sigma(given$sigma, K = length(theta_hat))
  labels <- .bef_site_labels(names(theta_hat), length(theta_hat))
  limits <- .bef_resolution_rule()
  L <- vapply(L, function(x) {
    .bef_assert_integer_scalar(x, "L", lower = limits$L_min, upper = limits$L_max)
  }, integer(1))
  M <- vapply(M, function(x) .bef_assert_integer_scalar(x, "M", lower = 3L, upper = 10L),
              integer(1))
  expansion <- vapply(expansion, function(x) .bef_assert_number(x, "expansion", lower = 0, upper = 5),
                      numeric(1))
  grid_method <- vapply(grid_method, .bef_validate_grid_method, character(1), USE.NAMES = FALSE)
  component <- vapply(component, function(x) {
    .bef_match_grid_choice(x, c("point", "normal"), "component")
  }, character(1), USE.NAMES = FALSE)
  supports <- if (is.null(support)) {
    NULL
  } else {
    lapply(if (is.list(support)) support else list(support), .bef_assert_support)
  }
  if (!is.logical(keep_fits) || length(keep_fits) != 1L || is.na(keep_fits)) {
    .bef_abort_invalid_args("`keep_fits` must be TRUE or FALSE.", arg = "keep_fits")
  }
  sites <- .bef_comparison_sites(sites, theta_hat, labels)
  dots <- list(...)

  specs <- if (is.null(supports)) {
    expand.grid(L = L, M = M, expansion = expansion, grid_method = grid_method,
                component = component, support = NA_integer_,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    expand.grid(L = L, M = M, expansion = NA_real_, grid_method = "support",
                component = component, support = seq_along(supports),
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  n <- nrow(specs)
  site_columns <- paste0("site_", labels[sites])

  rows <- vector("list", n)
  loos <- vector("list", n)
  g_list <- vector("list", n)
  site_rows <- vector("list", n)
  fits <- vector("list", n)
  for (i in seq_len(n)) {
    spec <- specs[i, ]
    described <- .bef_comparison_describe(spec, supports)
    rlang::inform(sprintf("Fit %d of %d: %s", i, n, described), class = "bef_progress")
    args <- c(
      list(theta_hat = theta_hat, sigma = sigma, L = spec$L, M = spec$M,
           component = spec$component),
      if (is.na(spec$support)) {
        list(grid_method = spec$grid_method, expansion = spec$expansion)
      } else {
        list(support = supports[[spec$support]])
      },
      dots
    )
    fit <- tryCatch(
      suppressWarnings(do.call(bayes_efron_fit, args)),
      error = function(err) {
        .bef_warn(
          sprintf("The fit failed for %s: %s", described, conditionMessage(err)),
          "bef_grid_comparison_fit_failed",
          spec = i
        )
        NULL
      }
    )
    row <- data.frame(
      spec = i, L = spec$L, M = spec$M, expansion = spec$expansion,
      grid_method = spec$grid_method, component = spec$component,
      lower = NA_real_, upper = NA_real_,
      elpd_loo = NA_real_, se_elpd_loo = NA_real_, p_loo = NA_real_,
      mean_g = NA_real_, sd_g = NA_real_,
      stringsAsFactors = FALSE
    )
    row[site_columns] <- NA_real_
    row <- cbind(row, data.frame(
      max_rhat = NA_real_, min_ess_bulk = NA_real_, divergences = NA_real_,
      max_pareto_k = NA_real_, end_mass = NA_real_, adequate = NA,
      sampler_seconds = NA_real_
    ))
    if (!is.null(fit)) {
      lo <- suppressWarnings(loo::loo(fit))
      d <- diagnose(fit)
      table <- site_summary(fit, sites = sites)
      loos[[i]] <- lo
      g_list[[i]] <- g_summary(fit, type = "density")
      site_rows[[i]] <- cbind(spec = i, table)
      if (isTRUE(keep_fits)) {
        fits[[i]] <- fit
      }
      grid <- fit$metadata$data_list$grid
      record <- attr(fit$metadata, "grid", exact = TRUE)
      row$lower <- grid[1L]
      row$upper <- grid[length(grid)]
      row$elpd_loo <- lo$estimates["elpd_loo", "Estimate"]
      row$se_elpd_loo <- lo$estimates["elpd_loo", "SE"]
      row$p_loo <- lo$estimates["p_loo", "Estimate"]
      row$mean_g <- mean(fit$posterior$mean_g)
      row$sd_g <- mean(fit$posterior$sd_g)
      row[site_columns] <- as.list(table$mean)
      row$max_rhat <- suppressWarnings(max(d$rhat, na.rm = TRUE))
      row$min_ess_bulk <- suppressWarnings(min(d$ess_bulk, na.rm = TRUE))
      row$divergences <- sum(d$divergences, na.rm = TRUE)
      row$max_pareto_k <- suppressWarnings(max(lo$diagnostics$pareto_k, na.rm = TRUE))
      row$end_mass <- if (is.null(d$grid)) NA_real_ else max(d$grid$ends$mass)
      row$adequate <- isTRUE(record$resolution$adequate)
      row$sampler_seconds <- fit$metadata$runtime_seconds
    }
    rows[[i]] <- row
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
  out <- out[, c("spec", "L", "M", "expansion", "grid_method", "component", "lower", "upper",
                 "elpd_loo", "se_elpd_loo", "p_loo", "elpd_diff", "se_diff",
                 "mean_g", "sd_g", site_columns,
                 "max_rhat", "min_ess_bulk", "divergences", "max_pareto_k",
                 "end_mass", "adequate", "sampler_seconds")]
  names(g_list) <- as.character(seq_len(n))
  structure(
    out,
    class = c("bef_grid_comparison", "data.frame"),
    g = g_list,
    sites = do.call(rbind, site_rows),
    fits = if (isTRUE(keep_fits)) stats::setNames(fits, as.character(seq_len(n))) else NULL
  )
}

# The sites whose summaries a comparison keeps, as positions: those asked
# for by position or by label, or the two with the extreme estimates.
.bef_comparison_sites <- function(sites, theta_hat, labels) {
  K <- length(theta_hat)
  if (is.null(sites)) {
    return(unique(c(which.min(theta_hat), which.max(theta_hat))))
  }
  if (is.character(sites)) {
    at <- match(sites, labels)
    if (anyNA(at)) {
      .bef_abort_invalid_args(
        sprintf("`sites` names no site: %s.",
                paste(sprintf("\"%s\"", sites[is.na(at)]), collapse = ", ")),
        arg = "sites"
      )
    }
    return(unique(at))
  }
  if (!is.numeric(sites) || length(sites) == 0L || anyNA(sites) ||
      any(sites != floor(sites)) || any(sites < 1L) || any(sites > K)) {
    .bef_abort_invalid_args(
      sprintf("`sites` must hold positions between 1 and %d, or labels of sites.", K),
      arg = "sites"
    )
  }
  unique(as.integer(sites))
}

# One specification in words, for the progress line and for a warning.
.bef_comparison_describe <- function(spec, supports) {
  range <- if (is.na(spec$support)) {
    sprintf("rule \"%s\", expansion %s", spec$grid_method, format(spec$expansion))
  } else {
    sprintf("support %s to %s", format(signif(supports[[spec$support]][1L], 4)),
            format(signif(supports[[spec$support]][2L], 4)))
  }
  sprintf("L = %d, M = %d, %s%s", spec$L, spec$M, range,
          if (identical(spec$component, "normal")) ", normal components" else "")
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
  attr(show, "g") <- attr(show, "sites") <- attr(show, "fits") <- NULL
  # Columns that say the same for every row are named once above the table.
  constant <- intersect(c("grid_method", "component", "expansion", "L", "M"), names(show))
  constant <- constant[vapply(constant, function(name) {
    nrow(show) > 1L && length(unique(show[[name]])) == 1L && !anyNA(show[[name]])
  }, logical(1))]
  if (length(constant) > 0L) {
    cat("In every row: ", paste(sprintf("%s = %s", constant, vapply(constant, function(name) {
      format(show[[name]][[1L]])
    }, character(1))), collapse = ", "), "\n\n", sep = "")
    show <- show[setdiff(names(show), constant)]
  }
  show <- show[setdiff(names(show), c("spec", "sampler_seconds", "min_ess_bulk", "p_loo"))]
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], function(v) signif(v, 3))
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
  if (!is.null(x$end_mass)) {
    leaning <- which(!is.na(x$end_mass) & x$end_mass > .bef_grid_ends_rule()$threshold)
    if (length(leaning) > 0L) {
      cat("\nRows in which a site has more than 5% of its posterior in the outer 2% of the grid:",
          paste(leaning, collapse = ", "), "\n")
    }
  }
  ok <- which(!is.na(x$elpd_diff) & !is.na(x$se_diff))
  if (length(ok) > 1L && all(abs(x$elpd_diff[ok]) <= 2 * x$se_diff[ok])) {
    cat("\nelpd_loo does not separate these specifications: all are within two standard",
        "\nerrors of the best. Compare sd_g and the site columns.\n")
  }
  invisible(x)
}
