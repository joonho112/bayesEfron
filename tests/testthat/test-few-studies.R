# The table of repeated samples with 19 studies that
# vignette("teacher-expectancy") reads: its form, and that it describes the
# stored fit of those studies.

few_studies <- function() {
  utils::read.csv(system.file("verification", "few-studies.csv", package = "bayesEfron"),
                  stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("the table has a row for each distribution and analysis", {
  table <- few_studies()
  scenarios <- c("normal, SD 0.14", "normal, SD 0.31", "fitted g")
  methods <- c("bayesEfron", "bayesEfron, adapt_delta = 0.99", "normal", "known g")
  expect_identical(table$scenario, rep(scenarios, each = 4L))
  expect_identical(table$method, rep(methods, times = 3L))
  expect_true(all(table$replications == 200L))
  shares <- c("sd_zero", "inclusion", "largest_inclusion", "own_as_large", "reml_as_small",
              "both_as_observed")
  for (column in shares) {
    expect_true(all(table[[column]] >= 0 & table[[column]] <= 1), info = column)
  }
  # The diagnostics of the sampler belong to the fits of this package.
  own <- grepl("^bayesEfron", table$method)
  expect_false(anyNA(table$with_divergences[own]))
  expect_true(all(is.na(table$with_divergences[!own])))
  expect_true(all(c("generated_on", "cmdstan_version") %in% names(table)))
})

test_that("the table starts from the stored fit of the 19 studies", {
  table <- few_studies()
  fit <- bayesEfron::raudenbush_fit
  sd_g <- diagnose(fit)$grid$sd_g[["mean"]]
  reml <- .bef_reml_normal(bayesEfron::raudenbush1985$yi, bayesEfron::raudenbush1985$sei)
  expect_equal(unique(table$observed_own_sd), sd_g, tolerance = 1e-10)
  # REML is an iteration, and its result on another platform agrees with the
  # stored one to about nine digits.
  expect_equal(unique(table$observed_reml_sd), reml$sd, tolerance = 1e-8)
  truth <- unique(table[, c("scenario", "true_sd")])
  expect_equal(truth$true_sd[1:2], c(reml$sd, sd_g), tolerance = 1e-8)
  # The third distribution is the posterior mean of g, a grid point moved
  # uniformly within its cell.
  mass <- g_summary(fit, type = "mass")
  probability <- mass$mean / sum(mass$mean)
  center <- sum(probability * mass$value)
  spacing <- diff(mass$value[1:2])
  expect_equal(truth$true_sd[3],
               sqrt(sum(probability * (mass$value - center)^2) + spacing^2 / 12),
               tolerance = 1e-10)
  # With the true distribution known, the estimate of its standard deviation
  # is the truth in every data set.
  known <- table[table$method == "known g", ]
  expect_equal(known$sd_mean, known$true_sd, tolerance = 1e-10)
})
