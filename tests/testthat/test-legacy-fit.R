legacy_ns <- function(name) getFromNamespace(name, "bayesEfron")

legacy_fit <- function(version) {
  readRDS(test_path("_fixtures", paste0("legacy-", version, ".rds")))
}

legacy_notice <- function(expr) {
  notices <- character()
  value <- withCallingHandlers(
    expr,
    bef_legacy_fit = function(w) {
      notices <<- c(notices, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, notices = notices)
}

# The summaries of the sites computed here from the draws of the saved fit,
# and not by the code of the package that computes them again.
legacy_reference <- function(fit) {
  dm <- as.matrix(posterior::as_draws_matrix(fit$draws))
  data <- fit$metadata$data_list
  site_matrix <- function(name) dm[, paste0(name, "[", seq_len(data$K), "]"), drop = FALSE]
  means <- site_matrix("theta_mean")
  sds <- site_matrix("theta_sd")
  g <- exp(dm[, paste0("log_g[", seq_len(data$L), "]"), drop = FALSE])
  modes <- vapply(seq_len(data$K), function(site) {
    probabilities <- sweep(g, 2L, stats::dnorm(
      data$theta_hat[site], data$grid, data$sigma[site]
    ), "*")
    probabilities <- probabilities / rowSums(probabilities)
    data$grid[which.max(colMeans(probabilities))]
  }, numeric(1))
  list(
    sd = unname(sqrt(colMeans(sds^2) + apply(means, 2L, stats::var))),
    mode = modes,
    effective = unname(rowSums(sweep(sds^2, 2L, data$sigma^2, "/"))),
    theta_rep = site_matrix("theta_rep")
  )
}

test_that("a fit saved by version 0.1 or 0.2 is summarized with the definitions of this version, and warns once per call", {
  for (version in c("v01", "v02")) {
    fit <- legacy_fit(version)
    unchanged <- fit
    reference <- legacy_reference(fit)
    # The table that this version computes from the draws of g in the fit.
    current <- legacy_ns(".bef_theta_summary")(legacy_ns(".bef_engine")(fit))
    result <- legacy_notice(as.data.frame(fit))
    expect_length(result$notices, 1L)
    expect_match(result$notices, "definitions of this version")
    expect_equal(result$value$sd, current$sd, tolerance = 1e-12)
    # The standard deviation of version 0.3.0 divides the variance between
    # draws by one less than their number.
    expect_equal(result$value$sd, reference$sd, tolerance = 0.02)
    expect_equal(result$value$mode, reference$mode, tolerance = 1e-12)
    expect_identical(result$value$map, result$value$mode)
    expect_equal(result$value$mean, fit$metadata$theta_summary$mean, tolerance = 1e-6)
    expect_equal(result$value$lower, current$lower)
    expect_equal(result$value$upper, current$upper)
    expect_identical(result$value$hpdi_lower, result$value$lower)

    result <- legacy_notice(summary(fit, level = 0.95))
    expect_length(result$notices, 1L)
    expect_equal(result$value$theta_summary$sd, current$sd)
    expect_equal(result$value$theta_summary$map, reference$mode)
    # At a level other than the stored one the limits come from the posterior
    # distribution; the quantiles of the theta_rep draws are still available
    # with method = "draws".
    exact <- suppressWarnings(confint(fit, level = 0.95))
    expect_equal(result$value$theta_summary$hpdi_lower, exact$lower)
    expect_equal(result$value$theta_summary$hpdi_upper, exact$upper)
    expected <- apply(reference$theta_rep, 2L, stats::quantile, c(0.025, 0.975))
    from_draws <- suppressWarnings(confint(fit, level = 0.95, method = "draws"))
    expect_equal(from_draws$lower, unname(expected[1L, ]))
    expect_equal(from_draws$upper, unname(expected[2L, ]))
    expect_identical(attr(result$value, "summary_definition_version"), 2L)
    expect_warning(format(result$value), NA)

    for (method in list(
      function() coef(fit, type = "map"),
      function() vcov(fit),
      function() confint(fit, level = 0.8),
      function() fitted(fit),
      function() format(fit, use_cli = FALSE),
      function() capture.output(print(fit, use_cli = FALSE))
    )) {
      expect_length(legacy_notice(method())$notices, 1L)
    }
    expect_equal(unname(legacy_notice(coef(fit, type = "map"))$value), reference$mode)
    expect_equal(unname(diag(legacy_notice(vcov(fit))$value)), current$sd^2)
    expect_identical(fit, unchanged)
  }
})

test_that("a fit saved by version 0.1 has its effective parameters computed again, and diagnose() and predict() refuse it", {
  fit <- legacy_fit("v01")
  reference <- legacy_reference(fit)
  repaired <- legacy_notice(legacy_ns(".bef_prepare_fit")(fit))$value
  effective <- repaired$posterior$effective_params
  expect_equal(effective, reference$effective, tolerance = 2e-6)
  expect_equal(effective, fit$metadata$data_list$K - fit$posterior$effective_params)
  dm <- as.matrix(posterior::as_draws_matrix(legacy_notice(posterior::as_draws(fit))$value))
  expect_equal(as.numeric(dm[, "effective_params"]), effective)
  expect_equal(legacy_notice(summary(fit))$value$diagnostics$effective_params$mean, mean(effective))
  expect_equal(attr(legacy_notice(logLik(fit))$value, "df"), mean(effective))
  expect_equal(as.numeric(legacy_notice(logLik(fit))$value),
               fit$metadata$log_marginal_likelihood_summary$mean)
  expect_true(legacy_ns(".bef_is_v01_fit")(repaired))
  expect_error(diagnose(repaired), "Refit", class = "bef_invalid_fit")
  expect_error(predict(repaired), "Refit", class = "bef_invalid_fit")
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_error(loo::loo(repaired), "Refit", class = "bef_invalid_fit")
    expect_error(loo::waic(repaired), "Refit", class = "bef_invalid_fit")
  }
})

test_that("plots of a fit saved by version 0.1 or 0.2 give one warning each with either backend, or ask for a refit", {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  for (version in c("v01", "v02")) {
    fit <- legacy_fit(version)
    for (type in c("caterpillar", "g", "diagnostic")) {
      result <- legacy_notice(plot(fit, type = type, backend = "base"))
      expect_length(result$notices, 1L)
      expect_null(result$value)
      if (identical(version, "v01") && identical(type, "diagnostic")) {
        # The ggplot2 plot draws what diagnose() reports, which a fit of
        # version 0.1 does not hold.
        expect_error(suppressWarnings(plot(fit, type = type)), "Refit",
                     class = "bef_invalid_fit")
      } else {
        result <- legacy_notice(plot(fit, type = type, backend = "ggplot2"))
        expect_length(result$notices, 1L)
        expect_s3_class(result$value, "ggplot")
      }
    }
  }
})

test_that("a fit made by this version keeps its summaries and gives no warning, even without its version attribute", {
  fit <- raudenbush_fit
  attr(fit$metadata, "summary_definition_version") <- NULL
  unchanged <- fit
  result <- legacy_notice(legacy_ns(".bef_prepare_fit")(fit))
  expect_length(result$notices, 0L)
  attr(result$value$metadata, "summary_definition_version") <- NULL
  expect_identical(result$value, fit)
  expect_warning(summary(fit), NA)
  expect_warning(as.data.frame(fit), NA)
  expect_identical(fit, unchanged)
})

test_that("the methods refuse a fit saved by version 0.2 that lacks its draws of g and a fit made by a later version", {
  incomplete <- legacy_fit("v02")
  keep <- !grepl("^(log_g|g)\\[", dimnames(incomplete$draws)[[3L]])
  incomplete$draws <- incomplete$draws[, , keep, drop = FALSE]
  future <- legacy_fit("v02")
  attr(future$metadata, "summary_definition_version") <- 3L
  expect_error(diagnose(future), "later version", class = "bef_invalid_fit")
  for (method in list(summary, coef, vcov, as.data.frame, logLik, format,
                      posterior::as_draws, function(x) confint(x),
                      function(x) plot(x, backend = "base"),
                      function(x) capture.output(print(x)))) {
    expect_error(method(incomplete), "Fit the model again", class = "bef_invalid_fit")
    expect_error(method(future), "later version", class = "bef_invalid_fit")
  }
})

test_that("a summary saved by this version prints the same, and one saved by an earlier or a later version is an error", {
  current <- summary(raudenbush_fit)
  serialized <- unserialize(serialize(current, NULL))
  expect_identical(format(serialized, use_cli = FALSE), format(current, use_cli = FALSE))
  attr(serialized, "summary_definition_version") <- NULL
  expect_error(format(serialized), "again on the fit", class = "bef_invalid_fit")
  expect_error(print(serialized), "again on the fit", class = "bef_invalid_fit")
  # A summary saved by version 0.3.0 has summary_definition_version 1.
  attr(serialized, "summary_definition_version") <- 1L
  expect_error(format(serialized), "again on the fit", class = "bef_invalid_fit")
  attr(serialized, "summary_definition_version") <- 3L
  expect_error(format(serialized), "later version", class = "bef_invalid_fit")
})
