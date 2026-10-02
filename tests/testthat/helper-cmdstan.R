# Skip a test unless CmdStan can compile and run a model on this machine.
skip_if_no_cmdstan <- function() {
  skip_on_cran()
  # The availability probe can warn when CMDSTAN names an unavailable path.
  # Keep a deliberate no-CmdStan skip quiet; warnings from actual fits are
  # checked by the tests themselves.
  suppressWarnings(skip_if_not_installed("cmdstanr"))

  path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(err) NA_character_
  )
  skip_if(
    length(path) != 1L || is.na(path) || !nzchar(path) || !dir.exists(path),
    "CmdStan is not installed."
  )

  toolchain <- tryCatch(
    suppressMessages(suppressWarnings(
      cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
    )),
    error = function(err) FALSE
  )
  skip_if(
    identical(toolchain, FALSE),
    "The C++ toolchain that CmdStan needs is not available."
  )

  invisible(TRUE)
}

# Share of sites whose true effect lies inside the equal-tailed posterior
# interval. `draws` has one column per site.
interval_coverage <- function(draws, theta_true, level = 0.9) {
  if (!is.matrix(draws)) {
    draws <- as.matrix(draws)
  }
  probs <- getFromNamespace(".bef_interval_probs", "bayesEfron")(level)
  intervals <- vapply(
    seq_len(ncol(draws)),
    function(site) {
      posterior::quantile2(draws[, site], probs = probs, names = FALSE)
    },
    numeric(2L)
  )
  mean(theta_true >= intervals[1L, ] & theta_true <= intervals[2L, ])
}
