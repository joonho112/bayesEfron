tier0_ns <- function(name) {
  getFromNamespace(name, "bayesEfron")
}

tier0_target <- function(target_id, statuses = "active") {
  .bef_target(target_id, statuses = statuses)
}

tier0_target_any_status <- function(target_id) {
  tier0_target(target_id, statuses = c("active", "deferred"))
}

tier0_object_sha256 <- function(x) {
  digest::digest(serialize(x, NULL, version = 3), algo = "sha256")
}

tier0_expect_current_runtime_env <- function(target_id) {
  recorded <- .bef_byte_runtime_env(target_id)
  current <- as.list(.bef_current_runtime_env(target_id)[1L, , drop = FALSE])
  expect_identical(
    stats::setNames(as.character(unlist(current[names(recorded)])), names(recorded)),
    stats::setNames(as.character(unlist(recorded)), names(recorded))
  )
}

tier0_grid_inputs <- function(grid_method) {
  fixture <- readRDS(
    .bef_testthat_file(
      file.path("_fixtures", "grid", paste0(grid_method, "_inputs.rds")),
      mustWork = TRUE
    )
  )
  if (!is.list(fixture) || !is.list(fixture$data)) {
    stop("Tier 0 grid input fixture must be a list with a `data` field.")
  }
  fixture$data
}

tier0_summary <- function(mean = 0) {
  list(
    mean = mean,
    sd = 0.1,
    q5 = mean - 0.1,
    q50 = mean,
    q95 = mean + 0.1
  )
}

tier0_fit <- function(K = 5L, S = 6L, L = 21L, M = 4L) {
  theta_rep <- matrix(seq(-0.5, 0.5, length.out = S * K), nrow = S, ncol = K)
  theta_mean <- matrix(rep(seq(-0.2, 0.2, length.out = K), each = S), nrow = S)
  theta_sd <- matrix(rep(seq(0.1, 0.2, length.out = K), each = S), nrow = S)
  theta_map <- matrix(rep(seq(-0.25, 0.15, length.out = K), each = S), nrow = S)
  theta_summary <- data.frame(
    site = seq_len(K),
    mean = colMeans(theta_mean),
    sd = colMeans(theta_sd),
    hpdi_lower = apply(theta_rep, 2L, stats::quantile, probs = 0.05, names = FALSE),
    hpdi_upper = apply(theta_rep, 2L, stats::quantile, probs = 0.95, names = FALSE),
    map = colMeans(theta_map)
  )
  metadata <- list(
    model_family = "RE",
    grid_method = "paper_realdata",
    seed = 20260511L,
    cmdstan_version = "2.38.0",
    stan_file_sha256 = paste(rep("a", 64L), collapse = ""),
    data_list = list(
      K = as.integer(K),
      theta_hat = seq(-0.4, 0.4, length.out = K),
      sigma = seq(0.1, 0.3, length.out = K),
      L = as.integer(L),
      grid = seq(-1, 1, length.out = L),
      M = as.integer(M),
      B = matrix(seq_len(L * M) / 100, nrow = L, ncol = M)
    ),
    runtime_seconds = 1.5,
    mean_g_summary = tier0_summary(0),
    var_g_summary = tier0_summary(1),
    theta_summary = theta_summary,
    theta_rep_draws = theta_rep,
    effective_params_summary = tier0_summary(3),
    log_marginal_likelihood_summary = tier0_summary(-10)
  )
  attr(metadata, "sd_g_summary") <- tier0_summary(1)
  attr(metadata, "diagnostics") <- list(
    rhat = 1.01,
    ess_bulk = 500,
    ess_tail = 450,
    divergences = 0,
    max_treedepth = 0
  )
  attr(metadata, "diagnostic_skipped") <- character()
  attr(metadata, "sampler_diagnostics_failed") <- character()

  variables <- c("mean_g", "var_g", "sd_g", paste0("g[", seq_len(L), "]"))
  tier0_ns("validate_bef_fit_re")(
    tier0_ns("new_bef_fit_re")(
      draws = array(
        seq_len(S * length(variables)) / 100,
        dim = c(S, 1L, length(variables)),
        dimnames = list(NULL, NULL, variables)
      ),
      metadata = metadata,
      posterior = list(
        mean_g = seq(-0.1, 0.1, length.out = S),
        var_g = seq(0.8, 1.2, length.out = S),
        sd_g = seq(0.9, 1.1, length.out = S),
        theta_map = theta_map,
        theta_mean = theta_mean,
        theta_sd = theta_sd,
        theta_rep = theta_rep,
        effective_params = rep(3, S),
        log_marginal_likelihood = rep(-10, S)
      )
    )
  )
}

test_that("Tier 0 softmax target preserves simplex normalization", {
  target <- tier0_target("tier0_softmax_sum_to_one")
  grid <- make_efron_grid(c(-1, 0, 2, 3), c(0.2, 0.3, 0.4, 0.5), M = 4L)
  alpha <- c(-0.25, 0.1, 0.35, -0.05)

  log_w <- drop(grid$B %*% alpha)
  shifted <- log_w - max(log_w)
  g <- exp(shifted) / sum(exp(shifted))

  expect_equal(sum(g), as.numeric(target$expected_value), tolerance = target$tolerance_value)
  expect_true(all(is.finite(g)))
  expect_true(all(g >= 0))

  stan_source <- readLines(.bef_source_file("inst", "stan", "efron_re.stan"), warn = FALSE)
  expect_true(any(grepl("log_softmax\\(log_w\\)", stan_source)))
  expect_true(any(grepl("simplex\\[L\\] g = softmax\\(log_w\\)", stan_source)))
})

test_that("Tier 0 grid, Stan-data, metadata, and GQ contracts are ledger-bound", {
  augmented <- tier0_target("tier0_augmented_rank")
  no_nan <- tier0_target("tier0_b_matrix_no_nan")
  k_target <- tier0_target("tier0_prepare_stan_data_K_invariance")
  metadata_target <- tier0_target("tier0_metadata_field_count")
  gq_target <- tier0_target("tier0_gq_count")

  theta_hat <- c(-1, 0, 2, 3)
  sigma <- c(0.2, 0.3, 0.4, 0.5)
  grid <- make_efron_grid(theta_hat, sigma, M = 4L)
  stan_data <- tier0_ns("prepare_stan_data")(list(theta_hat = theta_hat, sigma = sigma), grid)

  expect_true(all(is.finite(grid$B)), info = no_nan$target_id)
  expect_equal(
    qr(cbind(1, grid$B), tol = sqrt(.Machine$double.eps))$rank,
    grid$M + 1L,
    tolerance = augmented$tolerance_value,
    info = augmented$target_id
  )
  expect_equal(stan_data$K, length(theta_hat), tolerance = k_target$tolerance_value)
  expect_equal(length(tier0_ns(".bef_fit_re_metadata_fields")()), 13L)
  expect_equal(
    length(tier0_fit()$metadata),
    as.integer(metadata_target$expected_value),
    tolerance = metadata_target$tolerance_value
  )
  expect_equal(
    length(tier0_ns(".bef_generated_quantity_fields")()),
    as.integer(gq_target$expected_value),
    tolerance = gq_target$tolerance_value
  )
})

test_that("Tier 0 posterior finiteness contracts are ledger-bound", {
  theta_sd_finite <- tier0_target("tier0_theta_sd_finite")
  theta_sd_nonneg <- tier0_target("tier0_theta_sd_nonneg")
  effective_params_finite <- tier0_target("tier0_effective_params_finite")
  effective_params_diagnostic <- tier0_target("tier0_effective_params_in_K_range_diagnostic")
  fit <- tier0_fit()

  expect_true(all(is.finite(fit$posterior$theta_sd)), info = theta_sd_finite$target_id)
  expect_true(all(fit$posterior$theta_sd >= 0), info = theta_sd_nonneg$target_id)
  expect_true(
    all(is.finite(fit$posterior$effective_params)),
    info = effective_params_finite$target_id
  )
  expect_true(
    all(is.finite(fit$metadata$effective_params_summary$mean)),
    info = effective_params_finite$target_id
  )
  expect_false(
    isTRUE(effective_params_diagnostic$release_blocking),
    info = effective_params_diagnostic$target_id
  )
  expect_true(
    fit$metadata$effective_params_summary$mean >= 0 &&
      fit$metadata$effective_params_summary$mean <= fit$metadata$data_list$K,
    info = effective_params_diagnostic$target_id
  )
})

test_that("Tier 0 paper_realdata grid BYTE fixture is reproducible", {
  target_id <- "tier0_grid_method_realdata_byte_identity"
  target <- tier0_target_any_status(target_id)
  expect_identical(target$tolerance_type, "BYTE")
  tier0_expect_current_runtime_env(target_id)

  expected <- readRDS(.bef_target_fixture_path(target_id, mustWork = TRUE))
  inputs <- tier0_grid_inputs("paper_realdata")
  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    expansion = expected$expansion,
    grid_method = "paper_realdata"
  )

  expect_identical(actual, expected)
  expect_identical(tier0_object_sha256(actual), tier0_object_sha256(expected))

  wrong_expansion <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    expansion = 0.25,
    grid_method = "paper_realdata"
  )
  expect_false(identical(wrong_expansion, expected))
  expect_false(identical(
    tier0_object_sha256(wrong_expansion),
    tier0_object_sha256(expected)
  ))
})

test_that("Tier 0 paper_simulation grid BYTE fixture is reproducible", {
  target_id <- "tier0_grid_method_simulation_byte_identity"
  target <- tier0_target_any_status(target_id)
  expect_identical(target$tolerance_type, "BYTE")
  tier0_expect_current_runtime_env(target_id)

  expected <- readRDS(.bef_target_fixture_path(target_id, mustWork = TRUE))
  inputs <- tier0_grid_inputs("paper_simulation")
  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_true
  )

  expect_identical(actual, expected)
  expect_identical(tier0_object_sha256(actual), tier0_object_sha256(expected))

  missing_oracle <- tryCatch(
    make_efron_grid(
      theta_hat = inputs$theta_hat,
      sigma = inputs$sigma,
      L = expected$L,
      M = expected$M,
      grid_method = "paper_simulation"
    ),
    error = identity
  )
  expect_s3_class(missing_oracle, "bef_err_grid_oracle_required")

  wrong_oracle <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_hat
  )
  expect_false(identical(wrong_oracle, expected))
  expect_false(identical(
    tier0_object_sha256(wrong_oracle),
    tier0_object_sha256(expected)
  ))
})

test_that("Tier 0 paper_sensitivity grid BYTE fixture is reproducible", {
  target_id <- "tier0_grid_method_sensitivity_byte_identity"
  target <- tier0_target_any_status(target_id)
  expect_identical(target$tolerance_type, "BYTE")
  tier0_expect_current_runtime_env(target_id)

  expected <- readRDS(.bef_target_fixture_path(target_id, mustWork = TRUE))
  inputs <- tier0_grid_inputs("paper_sensitivity")
  expect_identical(expected$expansion, inputs$bound_expansion)
  expect_identical(inputs$bound_expansion, 0.5)
  actual <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true,
    bound_expansion = inputs$bound_expansion
  )

  expect_identical(actual, expected)
  expect_identical(tier0_object_sha256(actual), tier0_object_sha256(expected))

  default_bound <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true
  )
  expect_identical(default_bound, expected)

  missing_oracle <- tryCatch(
    make_efron_grid(
      theta_hat = inputs$theta_hat,
      sigma = inputs$sigma,
      L = expected$L,
      M = expected$M,
      grid_method = "paper_sensitivity",
      bound_expansion = expected$expansion
    ),
    error = identity
  )
  expect_s3_class(missing_oracle, "bef_err_grid_oracle_required")

  wrong_bound <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_sensitivity",
    theta_true = inputs$theta_true,
    bound_expansion = 0.25
  )
  expect_false(identical(wrong_bound, expected))
  expect_false(identical(
    tier0_object_sha256(wrong_bound),
    tier0_object_sha256(expected)
  ))

  absolute_padding <- make_efron_grid(
    theta_hat = inputs$theta_hat,
    sigma = inputs$sigma,
    L = expected$L,
    M = expected$M,
    grid_method = "paper_simulation",
    theta_true = inputs$theta_true
  )
  expect_false(identical(absolute_padding, expected))
  expect_false(identical(
    tier0_object_sha256(absolute_padding),
    tier0_object_sha256(expected)
  ))
})

test_that("Tier 0 locked checksum targets are represented as BYTE runtime rows", {
  for (target_id in c(
    "tier0_locked_stan_sha256",
    "tier0_locked_grid_R_sha256",
    "tier0_locked_data_prep_R_sha256"
  )) {
    target <- tier0_target(target_id)
    runtime_env <- .bef_byte_runtime_env(target_id)
    expect_identical(target$tolerance_type, "BYTE")
    expect_identical(runtime_env$target_id, target_id)
  }
})
