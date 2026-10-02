# Summaries of a fit that an earlier version of the package saved.
.bef_prepare_fit <- function(x) {
  .bef_require_inherits(x, "bef_fit", "object", "bef_invalid_fit")
  version <- .bef_fit_summary_version(x)
  if (identical(version, 1L)) {
    return(x)
  }

  old_v01 <- .bef_is_v01_fit(x)
  repaired <- tryCatch(
    .bef_recompute_saved_summary(x, old_v01 = old_v01),
    error = function(err) {
      .bef_abort_invalid_fit(
        paste0(
          "This fit was saved by an earlier version of bayesEfron, and its ",
          "summaries cannot be brought up to date because the draws or the ",
          "data stored in it are incomplete. Fit the model again with ",
          "`bayes_efron_fit()`."
        ),
        parent = err
      )
    }
  )
  if (length(repaired$changed) > 0L) {
    rlang::warn(
      paste0(
        "This fit was saved by an earlier version of bayesEfron. The ",
        .bef_enumerate(repaired$changed),
        " shown were computed again from its draws, with the definitions of ",
        "this version. The saved object is unchanged."
      ),
      class = c("bef_legacy_fit", "bef_warning"),
      call = .bef_user_frame()
    )
  }
  repaired$fit
}

.bef_fit_summary_version <- function(x) {
  version <- attr(x$metadata, "summary_definition_version", exact = TRUE)
  if (!is.null(version) && !identical(version, 1L)) {
    .bef_abort_invalid_fit(
      paste0(
        "This fit was made by a later version of bayesEfron than the one ",
        "installed. Update the package, or fit the model again."
      )
    )
  }
  version
}

.bef_recompute_saved_summary <- function(x, old_v01) {
  data <- x$metadata$data_list
  K <- data$K
  L <- data$L
  if (!.bef_is_whole_number(K) || K < 1L ||
      !.bef_is_whole_number(L) || L < 1L ||
      !is.numeric(data$grid) || length(data$grid) != L ||
      any(!is.finite(data$grid)) || any(diff(data$grid) <= 0) ||
      !is.numeric(data$theta_hat) || length(data$theta_hat) != K ||
      any(!is.finite(data$theta_hat)) ||
      !is.numeric(data$sigma) || length(data$sigma) != K ||
      any(!is.finite(data$sigma)) || any(data$sigma <= 0)) {
    stop("The saved grid, estimates or standard errors are incomplete.")
  }
  .bef_validate_theta_summary(x$metadata$theta_summary, K, "bef_invalid_fit")
  dm <- .bef_draws_matrix(x$draws)
  theta_mean <- .bef_draws_vector(dm, "theta_mean", K)
  theta_sd <- .bef_draws_vector(dm, "theta_sd", K)
  if (any(!is.finite(theta_mean)) || any(!is.finite(theta_sd)) ||
      any(theta_sd < 0)) {
    stop("The saved conditional means or standard deviations are invalid.")
  }
  log_names <- paste0("log_g[", seq_len(L), "]")
  if (all(log_names %in% colnames(dm))) {
    log_g <- dm[, log_names, drop = FALSE]
  } else {
    g <- .bef_draws_vector(dm, "g", L)
    if (any(!is.finite(g)) || any(g < 0) || any(rowSums(g) <= 0)) {
      stop("The saved grid probabilities are invalid.")
    }
    log_g <- log(g)
  }
  if (anyNA(log_g) || any(log_g == Inf) || any(rowSums(is.finite(log_g)) == 0)) {
    stop("The saved log probabilities are invalid.")
  }

  sd <- .bef_posterior_sd(theta_mean, theta_sd)
  mode <- .bef_posterior_mode(log_g, data$theta_hat, data$sigma, data$grid)
  changed <- character()
  old_sd <- x$metadata$theta_summary$sd
  if (any(abs(old_sd - sd) > 1e-10 * pmax(1, abs(sd)))) {
    x$metadata$theta_summary$sd <- sd
    changed <- c(changed, "posterior standard deviations")
  }
  if (!identical(as.numeric(x$metadata$theta_summary$map), as.numeric(mode))) {
    x$metadata$theta_summary$map <- as.numeric(mode)
    changed <- c(changed, "posterior modes")
  }

  if (old_v01) {
    # Version 0.1 stored K minus the sum of conditional variance ratios.
    old_effective <- .bef_draws_scalar(dm, "effective_params")
    if (any(!is.finite(old_effective))) {
      stop("The saved effective-parameter draws are invalid.")
    }
    effective <- K - old_effective
    if (!identical(as.numeric(effective), as.numeric(old_effective))) {
      changed <- c(changed, "effective parameters")
    }
    x$draws[, , "effective_params"] <- array(
      effective, dim = dim(x$draws)[1:2]
    )
    x$posterior$effective_params <- effective
    x$metadata$effective_params_summary <- .bef_summary_vector(effective)
  }
  attr(x$metadata, "summary_definition_version") <- 1L
  list(fit = x, changed = changed)
}

# "a", "a and b", "a, b and c"
.bef_enumerate <- function(x) {
  if (length(x) < 2L) {
    return(x)
  }
  paste(paste(x[-length(x)], collapse = ", "), "and", x[[length(x)]])
}

# A saved summary holds no draws from which it could be computed again.
.bef_check_saved_summary <- function(x) {
  version <- attr(x, "summary_definition_version", exact = TRUE)
  if (identical(version, 1L)) {
    return(invisible(x))
  }
  .bef_abort_invalid_fit(
    if (is.null(version)) {
      paste0(
        "This summary was saved by an earlier version of bayesEfron and ",
        "cannot be brought up to date, because it holds no draws. Call ",
        "`summary()` again on the fit."
      )
    } else {
      paste0(
        "This summary was made by a later version of bayesEfron than the one ",
        "installed. Update the package, or call `summary()` again on the fit."
      )
    }
  )
}
