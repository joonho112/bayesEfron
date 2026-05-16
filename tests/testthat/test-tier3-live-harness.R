test_that("Tier 3 live env-var envelope is explicit", {
  withr::local_envvar(c(
    BAYESEFRON_RUN_LIVE = NA_character_,
    BAYESEFRON_RUN_FULL_LIVE = NA_character_,
    BAYESEFRON_TIER3_FULL_LIVE_MATRIX = NA_character_,
    BAYESEFRON_TIER3_OK_TO_REFIT = NA_character_
  ))
  expect_false(live_cmdstan_run_requested())
  expect_false(live_cmdstan_smoke_run_requested())
  expect_false(live_cmdstan_full_run_requested())
  expect_false(tier3_full_live_matrix_requested())
  expect_false(tier3_full_live_refit_allowed())
  expect_identical(tier3_live_k_values(), integer())
  expect_identical(tier3_sampler_seed(), 20260509L)

  withr::local_envvar(c(BAYESEFRON_RUN_LIVE = "1"))
  expect_true(live_cmdstan_run_requested())
  expect_true(live_cmdstan_smoke_run_requested())
  expect_false(live_cmdstan_full_run_requested())
  expect_identical(tier3_live_k_values(), 50L)

  withr::local_envvar(c(BAYESEFRON_RUN_FULL_LIVE = "1"))
  expect_true(live_cmdstan_run_requested())
  expect_true(live_cmdstan_smoke_run_requested())
  expect_true(live_cmdstan_full_run_requested())
  expect_identical(tier3_live_k_values(), tier3_release_k_values())

  withr::local_envvar(c(BAYESEFRON_TIER3_OK_TO_REFIT = "0"))
  expect_false(tier3_full_live_refit_allowed())
  withr::local_envvar(c(BAYESEFRON_TIER3_OK_TO_REFIT = "true"))
  expect_false(tier3_full_live_refit_allowed())
  withr::local_envvar(c(BAYESEFRON_TIER3_OK_TO_REFIT = "1"))
  expect_true(tier3_full_live_refit_allowed())

  replay_matrix <- tempfile(fileext = ".csv")
  writeLines("placeholder", replay_matrix, useBytes = TRUE)
  withr::defer(unlink(replay_matrix))
  withr::local_envvar(c(
    BAYESEFRON_TIER3_FULL_LIVE_MATRIX = replay_matrix,
    BAYESEFRON_TIER3_OK_TO_REFIT = NA_character_
  ))
  expect_true(tier3_full_live_matrix_requested())
  expect_false(tier3_full_live_refit_allowed())
})

test_that("Tier 3 release ledger rows are active after full-live acceptance", {
  targets <- .bef_load_targets()
  release_ids <- tier3_target_id(tier3_release_k_values())
  release_rows <- targets[match(release_ids, targets$target_id), , drop = FALSE]

  expect_equal(release_rows$status, rep("active", length(release_ids)))
  expect_true(all(release_rows$release_blocking))
  expect_equal(
    release_rows$fixture_path,
    sprintf("_fixtures/lee_sui_K%d.rds", tier3_release_k_values())
  )

  diagnostic <- .bef_target("tier3_theta_mean_coverage_diagnostic_K500", targets = targets)
  runtime <- .bef_target("tier3_paper_realdata_runtime_under_300s", targets = targets)
  expect_equal(diagnostic$status, "deferred")
  expect_false(diagnostic$release_blocking)
  expect_equal(runtime$status, "deferred")
  expect_false(runtime$release_blocking)
})

test_that("Tier 3 coverage helper computes equal-tailed theta_rep coverage", {
  draws <- matrix(
    c(
      -1.0, 0.1, 1.0,
      -0.8, 0.2, 1.2,
      -0.6, 0.3, 1.4,
      -0.4, 0.4, 1.6,
      -0.2, 0.5, 1.8
    ),
    nrow = 5L,
    ncol = 3L,
    byrow = TRUE
  )
  theta_true <- c(-0.7, 0.35, 3.0)

  expect_equal(tier3_coverage(draws, theta_true, level = 0.8), 2 / 3)
})

test_that("Tier 3 live replication selector keeps smoke cheap and full mode complete", {
  v1_fixture <- list(metadata = list(fixture_format_version = "v1"))
  v2_fixture <- tier3_load_fixture(50L)

  expect_length(tier3_live_replications(v1_fixture, full = FALSE), 1L)
  expect_length(tier3_live_replications(v2_fixture, full = FALSE), 1L)
  expect_identical(tier3_live_replications(v2_fixture, full = FALSE)[[1L]], v2_fixture)

  full_replications <- tier3_live_replications(v2_fixture, full = TRUE)
  expect_length(full_replications, 20L)
  expect_equal(vapply(full_replications, `[[`, integer(1), "replication"), seq_len(20L))

  expect_error(
    tier3_live_replications(v1_fixture, full = TRUE),
    "v2 fixture with all 20 replications"
  )
})

test_that("Tier 3 aggregate coverage helper preserves denominator accounting", {
  expect_equal(tier3_aggregate_coverage(c(0.8, 0.9, 1.0)), 0.9)
  expect_error(tier3_aggregate_coverage(numeric()), "non-empty numeric")
  expect_error(tier3_aggregate_coverage(c("0.8", "0.9")), "non-empty numeric")
  expect_error(tier3_aggregate_coverage(c(0.8, NA_real_)), "non-finite")
  expect_error(tier3_aggregate_coverage(c(0.8, 1.2)), "\\[0, 1\\]")

  records <- lapply(seq_len(20L), function(replication) {
    list(replication = replication, coverage = 0.9, n_sites = 50L, status = "ok", error = NA_character_)
  })
  summary <- tier3_replication_coverage_summary(records, expected_replications = 20L, expected_sites = 1000L)

  expect_equal(summary$coverage, 0.9)
  expect_equal(summary$replication_coverage, rep(0.9, 20L))
  expect_equal(summary$n_replications_expected, 20L)
  expect_equal(summary$n_replications_fit, 20L)
  expect_equal(summary$n_replications_failed, 0L)
  expect_equal(summary$n_sites, 1000L)

  failed_records <- records
  failed_records[[2L]]$status <- "failed"
  failed_records[[2L]]$error <- "fit failed"
  expect_error(
    tier3_replication_coverage_summary(failed_records, expected_replications = 20L, expected_sites = 1000L),
    "replication fits failed"
  )

  missing_status_records <- records
  missing_status_records[[2L]]$status <- NULL
  expect_error(
    tier3_replication_coverage_summary(missing_status_records, expected_replications = 20L, expected_sites = 1000L),
    "replication fits failed"
  )

  expect_error(
    tier3_replication_coverage_summary(records[-20L], expected_replications = 20L, expected_sites = 1000L),
    "Expected 20 Tier 3 replication coverage records"
  )

  duplicated_replication_records <- records
  duplicated_replication_records[[2L]]$replication <- 1L
  expect_error(
    tier3_replication_coverage_summary(duplicated_replication_records, expected_replications = 20L, expected_sites = 1000L),
    "ordered replications"
  )

  uneven_records <- records
  uneven_records[[1L]]$n_sites <- 0L
  uneven_records[[2L]]$n_sites <- 100L
  expect_error(
    tier3_replication_coverage_summary(uneven_records, expected_replications = 20L, expected_sites = 1000L),
    "per replication"
  )

  expect_error(
    tier3_replication_coverage_summary(records, expected_replications = 20L, expected_sites = 500L),
    "site-level coverage indicators"
  )
})

test_that("Tier 3 full-live evidence matrix validates release targets", {
  matrix <- data.frame(
    target_id = tier3_target_id(tier3_release_k_values()),
    K = tier3_release_k_values(),
    status = rep("active", length(tier3_release_k_values())),
    release_blocking = TRUE,
    fixture_path = sprintf("_fixtures/lee_sui_K%d.rds", tier3_release_k_values()),
    expected_lower = 0.87,
    expected_upper = 0.92,
    coverage = rep(0.89, length(tier3_release_k_values())),
    in_band = TRUE,
    n_replications_expected = 20L,
    n_replications_fit = 20L,
    n_replications_failed = 0L,
    n_sites = 20L * tier3_release_k_values(),
    evidence_path = normalizePath(testthat::test_path("."), mustWork = TRUE)
  )
  path <- tempfile(fileext = ".csv")
  utils::write.csv(matrix, path, row.names = FALSE)

  loaded <- tier3_load_full_live_matrix(path)
  expect_equal(loaded$target_id, tier3_target_id(tier3_release_k_values()))
  expect_equal(loaded$K, tier3_release_k_values())

  bad <- matrix
  bad$expected_lower[1L] <- 0.9
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "lower bound does not match ledger",
    fixed = TRUE
  )

  bad <- matrix
  bad$expected_upper[1L] <- 0.95
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "upper bound does not match ledger",
    fixed = TRUE
  )

  bad <- matrix
  bad$status[1L] <- "deferred"
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "status does not match active ledger",
    fixed = TRUE
  )

  bad <- matrix
  bad$fixture_path[1L] <- "_fixtures/lee_sui_K50_missing.rds"
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "fixture_path does not match ledger",
    fixed = TRUE
  )

  bad <- matrix
  bad$coverage[1L] <- 0.5
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "coverage is outside the active ledger bounds",
    fixed = TRUE
  )

  bad <- matrix
  bad$in_band[1L] <- FALSE
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "includes out-of-band coverage",
    fixed = TRUE
  )

  bad <- matrix
  bad$n_replications_fit[1L] <- 19L
  bad$n_replications_failed[1L] <- 1L
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "replication accounting is invalid",
    fixed = TRUE
  )

  bad <- matrix
  bad$n_sites[1L] <- bad$n_sites[1L] + 1L
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "site denominator is invalid",
    fixed = TRUE
  )

  bad <- matrix
  bad$evidence_path[1L] <- file.path(tempdir(), "missing-tier3-evidence-dir")
  bad_path <- tempfile(fileext = ".csv")
  utils::write.csv(bad, bad_path, row.names = FALSE)
  expect_error(
    tier3_load_full_live_matrix(bad_path),
    "evidence_path does not exist",
    fixed = TRUE
  )
})

test_that("Tier 3 fit args use canonical Stan seed, not fixture provenance seed", {
  fixture <- list(
    theta_hat = c(-0.5, 0.5),
    sigma = c(0.2, 0.2),
    L = 51L,
    M = 6L,
    grid_method = "paper_realdata",
    expansion = 0.5,
    seed = 2026L
  )
  config <- list(chains = 1L, iter_warmup = 150L, iter_sampling = 20L)

  args <- tier3_fit_args(fixture, config)

  expect_identical(args$seed, 20260509L)
  expect_identical(args$seed, tier3_sampler_seed())
  expect_false(identical(args$seed, fixture$seed))
  expect_identical(fixture$seed, 2026L)
})

test_that("Tier 3 live smoke/full harness is gated by active fixtures", {
  skip_if_not(
    live_cmdstan_run_requested(),
    "Tier 3 live coverage is gated by BAYESEFRON_RUN_LIVE=1 or BAYESEFRON_RUN_FULL_LIVE=1."
  )

  K_values <- tier3_live_k_values()
  skip_if(length(K_values) == 0L, "No Tier 3 live K values selected.")

  matrix <- NULL
  if (live_cmdstan_full_run_requested() && tier3_full_live_matrix_requested()) {
    matrix <- tier3_load_full_live_matrix()
  }
  if (live_cmdstan_full_run_requested() && is.null(matrix) && !tier3_full_live_refit_allowed()) {
    fail(tier3_full_live_refit_guard_message())
    return(invisible())
  }
  if (is.null(matrix)) {
    live_cmdstan_skip_if_unavailable()
  }

  inactive <- K_values[!vapply(K_values, tier3_target_active, logical(1))]
  if (length(inactive) > 0L) {
    msg <- sprintf(
      "Tier 3 live coverage targets are not active: %s",
      paste(tier3_target_id(inactive), collapse = ", ")
    )
    if (live_cmdstan_full_run_requested()) {
      fail(msg)
      return(invisible())
    } else {
      skip(msg)
    }
  }

  missing_fixtures <- K_values[
    !file.exists(vapply(K_values, tier3_fixture_path, character(1), mustWork = FALSE))
  ]
  if (length(missing_fixtures) > 0L) {
    missing_paths <- vapply(
      missing_fixtures,
      tier3_fixture_path,
      character(1),
      mustWork = FALSE
    )
    msg <- sprintf(
      "Tier 3 live coverage fixtures are missing: %s",
      paste(missing_paths, collapse = ", ")
    )
    if (live_cmdstan_full_run_requested()) {
      fail(msg)
      return(invisible())
    } else {
      skip(msg)
    }
  }

  with_live_cmdstan_cache_root({
    for (K in K_values) {
      target <- .bef_target(tier3_target_id(K), statuses = "active")
      fixture <- tier3_load_fixture(K)
      tier3_validate_fixture(fixture, K)
      expected_replications <- if (live_cmdstan_full_run_requested()) 20L else 1L
      expected_sites <- if (live_cmdstan_full_run_requested()) {
        fixture$metadata$coverage$denominator
      } else {
        as.integer(K)
      }

      if (!is.null(matrix)) {
        matrix_row <- tier3_full_live_matrix_row(matrix, K)
        coverage <- matrix_row$coverage
        summary <- list(
          n_replications_expected = matrix_row$n_replications_expected,
          n_replications_fit = matrix_row$n_replications_fit,
          n_replications_failed = matrix_row$n_replications_failed,
          n_sites = matrix_row$n_sites
        )
      } else {
        config <- tier3_fit_config()
        replications <- tier3_live_replications(fixture)
        records <- lapply(replications, tier3_fit_replication_coverage, config = config)
        summary <- tier3_replication_coverage_summary(
          records,
          expected_replications = expected_replications,
          expected_sites = expected_sites
        )
        coverage <- summary$coverage
      }

      expect_equal(summary$n_replications_expected, expected_replications)
      expect_equal(summary$n_replications_fit, expected_replications)
      expect_equal(summary$n_replications_failed, 0L)
      expect_equal(summary$n_sites, expected_sites)

      if (live_cmdstan_full_run_requested()) {
        expect_gte(coverage, target$expected_value_lower)
        expect_lte(coverage, target$expected_value_upper)
      }
    }
  })
})
