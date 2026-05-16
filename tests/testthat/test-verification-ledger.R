test_that("verification target ledger has stable schema", {
  targets <- .bef_load_targets(refresh = TRUE)

  expect_named(targets, .bef_targets_schema())
  expect_equal(nrow(targets), 37L)
  expect_equal(anyDuplicated(targets$target_id), 0L)
  expect_true(all(targets$tier %in% 0:3))
  expect_true(all(targets$tolerance_type %in% c(
    "absolute", "relative", "interval", "BYTE"
  )))
  expect_true(all(targets$status %in% c("active", "deferred", "obsolete")))
  expect_true(all(is.na(targets$tolerance_value[
    targets$tolerance_type %in% c("interval", "BYTE")
  ])))
})

test_that("target lookup returns one row with parsed scalar fields", {
  target <- .bef_target("tier0_metadata_field_count", statuses = "active")

  expect_identical(target$target_id, "tier0_metadata_field_count")
  expect_identical(target$tier, 0L)
  expect_identical(target$tolerance_type, "absolute")
  expect_identical(target$release_blocking, TRUE)
  expect_identical(as.integer(target$expected_value), 13L)
})

test_that("runtime env rows cover BYTE targets exactly", {
  targets <- .bef_load_targets()
  runtime_env <- .bef_load_runtime_env(targets = targets, refresh = TRUE)
  byte_ids <- targets$target_id[targets$tolerance_type == "BYTE"]

  expect_named(runtime_env, .bef_runtime_env_schema())
  expect_setequal(runtime_env$target_id, byte_ids)
  expect_equal(anyDuplicated(runtime_env$target_id), 0L)
  expect_true(all(nzchar(runtime_env$R_version_string)))
  expect_true(all(nzchar(runtime_env$splines_pkg_version)))
  expect_true(all(nzchar(runtime_env$RDS_serialization_version)))
  expect_true(all(nzchar(runtime_env$platform)))
  expect_true(all(nzchar(runtime_env$BLAS_LAPACK)))
})

test_that("runtime env lookup is restricted to BYTE targets", {
  env <- .bef_byte_runtime_env("tier0_locked_stan_sha256")
  expect_identical(env$target_id, "tier0_locked_stan_sha256")

  expect_error(
    .bef_byte_runtime_env("tier0_metadata_field_count"),
    "not a BYTE target",
    fixed = TRUE
  )
})

test_that("fixture and source path helpers use source-tree-safe lookup", {
  source_path <- .bef_source_file("inst", "locked-core-checksums.txt")
  expect_true(file.exists(source_path))

  tier3_fixture <- .bef_target_fixture_path(
    "tier3_theta_rep_coverage_paper_realdata_K50",
    mustWork = FALSE
  )
  expect_match(tier3_fixture, "_fixtures/lee_sui_K50[.]rds")
})

test_that("release-mode ledger has no deferred release-blocking targets", {
  release_mode <- identical(Sys.getenv("BAYESEFRON_ENFORCE_RELEASE_LEDGER"), "1") ||
    identical(Sys.getenv("BAYESEFRON_RUN_FULL_LIVE"), "1")
  skip_if_not(
    release_mode,
    "Release-blocking ledger invariant is enforced only in release mode."
  )

  targets <- .bef_load_targets(refresh = TRUE)
  deferred_blockers <- targets$release_blocking & targets$status == "deferred"
  expect_false(
    any(deferred_blockers),
    info = paste(
      "Deferred release-blocking targets:",
      paste(targets$target_id[deferred_blockers], collapse = ", ")
    )
  )
})
