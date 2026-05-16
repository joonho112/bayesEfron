live_cmdstan_skip_if_unavailable <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  cmdstan_path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(err) NA_character_
  )
  skip_if(
    length(cmdstan_path) != 1L || is.na(cmdstan_path) ||
      !nzchar(cmdstan_path) || !dir.exists(cmdstan_path),
    "CmdStan is not configured for this test environment."
  )

  cmdstan_version <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(err) NA_character_
  )
  skip_if(
    length(cmdstan_version) != 1L || is.na(cmdstan_version) ||
      !nzchar(cmdstan_version),
    "CmdStan version is unavailable for this test environment."
  )

  toolchain_check <- tryCatch(
    suppressMessages(suppressWarnings(
      cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
    )),
    error = function(err) FALSE
  )
  skip_if_not(
    !identical(toolchain_check, FALSE),
    "CmdStan C++ toolchain is unavailable for this test environment."
  )

  invisible(list(path = cmdstan_path, version = cmdstan_version))
}

live_cmdstan_run_requested <- function() {
  live_cmdstan_smoke_run_requested() || live_cmdstan_full_run_requested()
}

live_cmdstan_smoke_run_requested <- function() {
  identical(Sys.getenv("BAYESEFRON_RUN_LIVE"), "1") ||
    live_cmdstan_full_run_requested()
}

live_cmdstan_full_run_requested <- function() {
  identical(Sys.getenv("BAYESEFRON_RUN_FULL_LIVE"), "1")
}

tier3_release_k_values <- function() {
  c(50L, 100L, 200L, 500L, 1500L)
}

tier3_live_k_values <- function() {
  if (live_cmdstan_full_run_requested()) {
    return(tier3_release_k_values())
  }
  if (live_cmdstan_smoke_run_requested()) {
    return(50L)
  }
  integer()
}

tier3_sampler_seed <- function() {
  20260509L
}

tier3_target_id <- function(K) {
  sprintf("tier3_theta_rep_coverage_paper_realdata_K%d", as.integer(K))
}

tier3_full_live_matrix_path <- function() {
  path <- Sys.getenv("BAYESEFRON_TIER3_FULL_LIVE_MATRIX", unset = "")
  if (!nzchar(path)) {
    return(NULL)
  }
  normalizePath(path, mustWork = TRUE)
}

tier3_full_live_matrix_requested <- function() {
  !is.null(tier3_full_live_matrix_path())
}

tier3_full_live_refit_allowed <- function() {
  identical(Sys.getenv("BAYESEFRON_TIER3_OK_TO_REFIT"), "1")
}

tier3_full_live_refit_guard_message <- function() {
  paste(
    "Tier 3 full-live fresh refit requires BAYESEFRON_TIER3_OK_TO_REFIT=1",
    "or BAYESEFRON_TIER3_FULL_LIVE_MATRIX."
  )
}

tier3_fixture_path <- function(K, mustWork = TRUE) {
  .bef_target_fixture_path(tier3_target_id(K), mustWork = mustWork)
}

tier3_load_fixture <- function(K) {
  readRDS(tier3_fixture_path(K, mustWork = TRUE))
}

tier3_fit_args <- function(fixture, config = tier3_fit_config()) {
  list(
    theta_hat = fixture$theta_hat,
    sigma = fixture$sigma,
    L = fixture$L %||% 101L,
    M = fixture$M %||% 6L,
    grid_method = fixture$grid_method %||% "paper_realdata",
    expansion = fixture$expansion %||% 0.5,
    chains = config$chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    seed = tier3_sampler_seed(),
    keep_cmdstan_fit = FALSE
  )
}

tier3_live_replications <- function(fixture, full = live_cmdstan_full_run_requested()) {
  if (isTRUE(full)) {
    if (!identical(fixture$metadata$fixture_format_version, "v2") || is.null(fixture$replications)) {
      stop("Full Tier 3 live mode requires a v2 fixture with all 20 replications.", call. = FALSE)
    }
    return(fixture$replications)
  }

  if (!is.null(fixture$replications)) {
    return(list(fixture))
  }

  list(fixture)
}

tier3_fixture_controls <- function(K) {
  controls <- data.frame(
    K = tier3_release_k_values(),
    L = c(51L, 71L, 81L, 101L, 101L),
    dataset_id = c(1L, 21L, 41L, 61L, 81L),
    seed = c(2026L, 2046L, 2066L, 2086L, 2106L),
    left_count = c(17L, 33L, 67L, 167L, 500L),
    right_count = c(33L, 67L, 133L, 333L, 1000L)
  )
  controls$dataset_id_end <- controls$dataset_id + 19L
  controls$seed_end <- controls$seed + 19L
  controls[match(as.integer(K), controls$K), , drop = FALSE]
}

tier3_replication_controls <- function(K) {
  controls <- tier3_fixture_controls(K)
  data.frame(
    K = as.integer(K),
    replication = seq_len(20L),
    dataset_id = seq.int(controls$dataset_id, controls$dataset_id_end),
    seed = seq.int(controls$seed, controls$seed_end),
    L = controls$L,
    M = 6L,
    grid_method = "paper_realdata",
    expansion = 0.5,
    bound_expansion = 0.5
  )
}

tier3_sha256_pattern <- function() {
  "^[0-9a-f]{64}$"
}

tier3_required_source_artifacts <- function() {
  c(
    "data_generation_script",
    "model_script",
    "public_small_k_script",
    "original_small_k_script",
    "stan_model",
    "simulation_rds",
    "archived_small_k_result"
  )
}

tier3_validate_source_artifacts <- function(artifacts) {
  expect_type(artifacts, "list")
  expect_true(all(tier3_required_source_artifacts() %in% names(artifacts)))

  for (name in tier3_required_source_artifacts()) {
    artifact <- artifacts[[name]]
    expect_type(artifact, "list")
    expect_true(all(c("role", "path", "sha256", "relevant_lines") %in% names(artifact)))
    expect_type(artifact$role, "character")
    expect_length(artifact$role, 1L)
    expect_true(nzchar(artifact$role))
    expect_type(artifact$path, "character")
    expect_length(artifact$path, 1L)
    expect_true(grepl("^/", artifact$path))
    expect_type(artifact$sha256, "character")
    expect_length(artifact$sha256, 1L)
    expect_match(artifact$sha256, tier3_sha256_pattern())
  }

  invisible(TRUE)
}

tier3_validate_fixture_vectors <- function(x, K, controls) {
  expect_length(x$theta_hat, K)
  expect_length(x$sigma, K)
  expect_length(x$theta_true, K)
  expect_type(x$theta_hat, "double")
  expect_type(x$sigma, "double")
  expect_type(x$theta_true, "double")
  expect_true(all(is.finite(x$theta_hat)))
  expect_true(all(is.finite(x$sigma)))
  expect_true(all(x$sigma > 0))
  expect_true(all(is.finite(x$theta_true)))
  expect_true(is.finite(diff(range(x$theta_hat))))
  expect_gt(diff(range(x$theta_hat)), 0)

  expect_length(x$source_row_id, K)
  expect_type(x$source_row_id, "integer")
  expect_equal(anyDuplicated(x$source_row_id), 0L)
  expect_true(all(x$source_row_id >= 1L & x$source_row_id <= 1500L))
  if (identical(as.integer(K), 1500L)) {
    expect_equal(x$source_row_id, seq_len(1500L))
  }

  expect_length(x$tower, K)
  expect_type(x$tower, "character")
  expect_setequal(unique(x$tower), c("left_theta_true_lt_0", "right_theta_true_ge_0"))
  expect_equal(sum(x$tower == "left_theta_true_lt_0"), controls$left_count)
  expect_equal(sum(x$tower == "right_theta_true_ge_0"), controls$right_count)

  invisible(TRUE)
}

tier3_validate_fixture_replication <- function(replication, K, controls, expected) {
  required_fields <- c(
    "replication", "dataset_id", "seed", "K",
    "theta_hat", "sigma", "theta_true", "source_row_id", "tower",
    "L", "M", "grid_method", "expansion", "bound_expansion", "grid"
  )

  expect_type(replication, "list")
  expect_true(all(required_fields %in% names(replication)))
  expect_equal(replication$replication, expected$replication)
  expect_equal(replication$dataset_id, expected$dataset_id)
  expect_equal(replication$seed, expected$seed)
  expect_equal(replication$K, as.integer(K))

  tier3_validate_fixture_vectors(replication, K, controls)

  expect_equal(replication$L, controls$L)
  expect_equal(replication$M, 6L)
  expect_equal(replication$grid_method, "paper_realdata")
  expect_equal(replication$expansion, 0.5)
  expect_equal(replication$bound_expansion, 0.5)

  expect_type(replication$grid, "list")
  expect_length(replication$grid$observed_range, 2L)
  expect_length(replication$grid$grid_bounds, 2L)
  expect_length(replication$grid$grid_expansion, 1L)
  expect_true(all(is.finite(replication$grid$observed_range)))
  expect_true(all(is.finite(replication$grid$grid_bounds)))
  expect_true(all(is.finite(replication$grid$grid_expansion)))
  observed_range <- range(replication$theta_hat)
  grid_expansion <- 0.5 * diff(observed_range)
  expect_equal(replication$grid$observed_range, as.numeric(observed_range))
  expect_equal(replication$grid$grid_bounds, as.numeric(c(
    observed_range[1L] - grid_expansion,
    observed_range[2L] + grid_expansion
  )))
  expect_equal(replication$grid$grid_expansion, as.numeric(grid_expansion))

  invisible(TRUE)
}

tier3_replication_payload_hash <- function(replication) {
  digest::digest(
    list(
      theta_hat = replication$theta_hat,
      sigma = replication$sigma,
      theta_true = replication$theta_true,
      source_row_id = replication$source_row_id,
      tower = replication$tower
    ),
    algo = "sha256"
  )
}

tier3_expect_unique_replication_payloads <- function(replications) {
  hashes <- vapply(replications, tier3_replication_payload_hash, character(1))
  expect_equal(anyDuplicated(hashes), 0L)
  invisible(hashes)
}

tier3_validate_fixture <- function(fixture, K) {
  common_required_fields <- c(
    "metadata", "K", "replication", "dataset_id", "seed",
    "theta_hat", "sigma", "theta_true", "source_row_id", "tower",
    "L", "M", "grid_method", "expansion", "bound_expansion"
  )
  controls <- tier3_fixture_controls(K)

  expect_type(fixture, "list")
  expect_true(all(common_required_fields %in% names(fixture)))
  expect_equal(as.integer(K), controls$K)

  metadata <- fixture$metadata
  expect_type(metadata, "list")
  expect_true("fixture_format_version" %in% names(metadata))
  expect_true(metadata$fixture_format_version %in% c("v1", "v2"))
  version <- metadata$fixture_format_version
  expect_equal(metadata$generator, "tools/generate-tier3-lee-sui-fixtures.R")
  expect_true("source" %in% names(metadata))
  expect_type(metadata$source, "list")
  expect_equal(metadata$source$scenario$I, 0.7)
  expect_equal(metadata$source$scenario$R, 9)
  expected_canonical_input <- if (identical(as.integer(K), 1500L)) {
    "part01_full_panel_simulation_rule"
  } else {
    "archived_appendix_small_k_result"
  }
  expect_equal(metadata$source$canonical_input, expected_canonical_input)
  expect_equal(metadata$source$archived_control_input, "archived_appendix_small_k_result")
  tier3_validate_source_artifacts(metadata$source$artifacts)
  expect_length(metadata$source$source_row_ids, K)
  expect_equal(metadata$source$source_row_ids, fixture$source_row_id)
  expect_length(metadata$source$first_source_row_ids, min(8L, as.integer(K)))
  expect_equal(metadata$source$first_source_row_ids, utils::head(fixture$source_row_id, 8L))

  expect_equal(metadata$K_scenarios, tier3_release_k_values())
  expect_equal(metadata$n_replications, 20L)
  expect_equal(metadata$replication, 1L)
  expect_equal(metadata$dataset_id, controls$dataset_id)
  expect_equal(metadata$seed_base, 2025L)
  expect_equal(metadata$seed_range, c(2026L, 2125L))
  expected_sampling_strategy <- if (identical(as.integer(K), 1500L)) {
    "part01_seeded_full_panel_resimulation"
  } else {
    "theta_true_split_at_0_proportional_towers_unsorted"
  }
  expect_equal(metadata$sampling_strategy, expected_sampling_strategy)
  expect_equal(metadata$RNGkind, c("Mersenne-Twister", "Inversion", "Rejection"))
  expect_type(metadata$grid, "list")
  expect_equal(metadata$grid$L, controls$L)
  expect_equal(metadata$grid$M, 6L)
  expect_equal(metadata$grid$grid_method, "paper_realdata")
  expect_equal(metadata$grid$expansion, 0.5)
  expect_equal(metadata$grid$bound_expansion, 0.5)

  if (identical(version, "v2")) {
    expect_true("replications" %in% names(fixture))
    expect_type(fixture$replications, "list")
    expect_length(fixture$replications, 20L)
    expect_equal(metadata$K, as.integer(K))
    expect_equal(metadata$primary_replication, 1L)
    expect_equal(metadata$replications, seq_len(20L))
    expect_equal(metadata$dataset_id_range, c(controls$dataset_id, controls$dataset_id_end))
    expect_equal(metadata$replication_seed_range, c(controls$seed, controls$seed_end))
    expect_equal(metadata$grid$per_replication_observed_range, TRUE)
    expect_type(metadata$coverage, "list")
    expect_equal(metadata$coverage$nominal_level, 0.9)
    expect_equal(metadata$coverage$aggregation, "mean_of_replication_site_coverages")
    expect_equal(metadata$coverage$n_replications, 20L)
    expect_equal(metadata$coverage$denominator, 20L * as.integer(K))
    expect_equal(metadata$coverage$failure_policy, "do_not_silently_drop_failed_replications")
    expect_equal(metadata$coverage$replication_identity_policy, "unique_input_payload_hashes")
    expect_equal(metadata$coverage$sampler_seed_policy, "tier3_sampler_seed")

    expected <- tier3_replication_controls(K)
    for (i in seq_along(fixture$replications)) {
      tier3_validate_fixture_replication(
        fixture$replications[[i]],
        K = K,
        controls = controls,
        expected = expected[i, , drop = FALSE]
      )
    }
    tier3_expect_unique_replication_payloads(fixture$replications)
  }

  expect_equal(fixture$K, as.integer(K))
  expect_equal(fixture$replication, 1L)
  expect_equal(fixture$dataset_id, controls$dataset_id)
  expect_equal(fixture$seed, controls$seed)
  tier3_validate_fixture_vectors(fixture, K, controls)

  expect_equal(fixture$L, controls$L)
  expect_equal(fixture$M, 6L)
  expect_equal(fixture$grid_method, "paper_realdata")
  expect_equal(fixture$expansion, 0.5)
  expect_equal(fixture$bound_expansion, 0.5)

  if (identical(version, "v2")) {
    primary <- fixture$replications[[1L]]
    expect_equal(fixture$theta_hat, primary$theta_hat)
    expect_equal(fixture$sigma, primary$sigma)
    expect_equal(fixture$theta_true, primary$theta_true)
    expect_equal(fixture$source_row_id, primary$source_row_id)
    expect_equal(fixture$tower, primary$tower)
    expect_equal(fixture$seed, primary$seed)
    expect_equal(fixture$dataset_id, primary$dataset_id)
    expect_equal(fixture$L, primary$L)
    expect_equal(fixture$M, primary$M)
    expect_equal(fixture$grid_method, primary$grid_method)
    expect_equal(fixture$expansion, primary$expansion)
    expect_equal(fixture$bound_expansion, primary$bound_expansion)
  }

  invisible(TRUE)
}

tier3_coverage <- function(draws, theta_true, level = 0.9) {
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

tier3_fit_replication_coverage <- function(replication, config = tier3_fit_config()) {
  tryCatch(
    {
      fit <- do.call(bayes_efron_fit, tier3_fit_args(replication, config))
      coverage <- tier3_coverage(fit$metadata$theta_rep_draws, replication$theta_true)
      list(
        replication = replication$replication %||% NA_integer_,
        coverage = coverage,
        n_sites = length(replication$theta_true),
        status = "ok",
        error = NA_character_
      )
    },
    error = function(err) {
      list(
        replication = replication$replication %||% NA_integer_,
        coverage = NA_real_,
        n_sites = length(replication$theta_true),
        status = "failed",
        error = conditionMessage(err)
      )
    }
  )
}

tier3_aggregate_coverage <- function(replication_coverages) {
  if (!is.numeric(replication_coverages) || length(replication_coverages) == 0L) {
    stop("Tier 3 replication coverage must be a non-empty numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(replication_coverages))) {
    stop("Tier 3 replication coverage contains non-finite values.", call. = FALSE)
  }
  if (any(replication_coverages < 0 | replication_coverages > 1)) {
    stop("Tier 3 replication coverage values must lie in [0, 1].", call. = FALSE)
  }

  mean(replication_coverages)
}

tier3_replication_coverage_summary <- function(records, expected_replications, expected_sites) {
  if (length(records) != expected_replications) {
    stop(sprintf(
      "Expected %d Tier 3 replication coverage records, got %d.",
      as.integer(expected_replications),
      length(records)
    ), call. = FALSE)
  }

  replication_ids <- vapply(records, function(record) as.integer(record$replication %||% NA_integer_), integer(1))
  expected_replication_ids <- seq_len(as.integer(expected_replications))
  if (!identical(replication_ids, expected_replication_ids)) {
    stop(sprintf(
      "Tier 3 replication coverage records must be ordered replications %s.",
      paste(expected_replication_ids, collapse = ", ")
    ), call. = FALSE)
  }

  status <- vapply(records, function(record) record$status %||% "missing", character(1))
  failed <- !identical(status, rep("ok", length(status)))
  if (failed) {
    errors <- vapply(records, function(record) record$error %||% NA_character_, character(1))
    failed_msg <- paste(
      sprintf("replication %s: %s", replication_ids[status != "ok"], errors[status != "ok"]),
      collapse = "; "
    )
    stop(sprintf("Tier 3 replication fits failed: %s", failed_msg), call. = FALSE)
  }

  coverage <- vapply(records, function(record) record$coverage, numeric(1))
  n_sites <- vapply(records, function(record) as.integer(record$n_sites), integer(1))
  total_sites <- sum(n_sites)
  expected_sites_per_replication <- as.integer(expected_sites / expected_replications)
  if (!identical(expected_sites_per_replication * as.integer(expected_replications), as.integer(expected_sites))) {
    stop("Tier 3 expected site denominator must divide evenly across replications.", call. = FALSE)
  }
  if (!all(n_sites == expected_sites_per_replication)) {
    stop(sprintf(
      "Expected %d Tier 3 site-level coverage indicators per replication.",
      expected_sites_per_replication
    ), call. = FALSE)
  }
  if (!identical(total_sites, as.integer(expected_sites))) {
    stop(sprintf(
      "Expected %d Tier 3 site-level coverage indicators, got %d.",
      as.integer(expected_sites),
      total_sites
    ), call. = FALSE)
  }
  if (any(!is.finite(coverage))) {
    stop("Tier 3 replication coverage contains non-finite values.", call. = FALSE)
  }

  list(
    coverage = tier3_aggregate_coverage(coverage),
    replication_coverage = coverage,
    n_replications_expected = as.integer(expected_replications),
    n_replications_fit = length(records),
    n_replications_failed = 0L,
    n_sites = total_sites
  )
}

tier3_full_live_matrix_schema <- function() {
  c(
    "target_id",
    "K",
    "status",
    "release_blocking",
    "fixture_path",
    "expected_lower",
    "expected_upper",
    "coverage",
    "in_band",
    "n_replications_expected",
    "n_replications_fit",
    "n_replications_failed",
    "n_sites",
    "evidence_path"
  )
}

tier3_full_live_matrix_path_is_absolute <- function(path) {
  grepl("^([A-Za-z]:[/\\\\]|/|~)", path)
}

tier3_full_live_matrix_root <- function(path) {
  matrix_dir <- dirname(normalizePath(path, mustWork = TRUE))
  candidates <- unique(c(
    dirname(dirname(matrix_dir)),
    bef_test_source_root(required = c("DESCRIPTION", "NAMESPACE"))
  ))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "DESCRIPTION"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  normalizePath(bef_test_source_root(required = c("DESCRIPTION", "NAMESPACE")), mustWork = TRUE)
}

tier3_resolve_full_live_matrix_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) {
    return(path)
  }
  path <- path.expand(path)
  if (tier3_full_live_matrix_path_is_absolute(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(root, path), winslash = "/", mustWork = FALSE)
}

tier3_load_full_live_matrix <- function(path = tier3_full_live_matrix_path(),
                                        targets = .bef_load_targets()) {
  if (is.null(path)) {
    stop("BAYESEFRON_TIER3_FULL_LIVE_MATRIX is not set.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Tier 3 full-live matrix does not exist.", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  matrix_root <- tier3_full_live_matrix_root(path)
  matrix <- utils::read.csv(
    path,
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(tier3_full_live_matrix_schema() %in% names(matrix))) {
    stop("Tier 3 full-live matrix is missing required columns.", call. = FALSE)
  }
  path_columns <- intersect(
    c("evidence_path", "evidence_log", "source_aggregate", "source_per_replication"),
    names(matrix)
  )
  for (column in path_columns) {
    matrix[[column]] <- vapply(
      matrix[[column]],
      tier3_resolve_full_live_matrix_path,
      character(1),
      root = matrix_root
    )
  }

  matrix <- matrix[, tier3_full_live_matrix_schema(), drop = FALSE]
  matrix$K <- as.integer(matrix$K)
  matrix$release_blocking <- .bef_parse_logical(matrix$release_blocking, "release_blocking")
  matrix$in_band <- .bef_parse_logical(matrix$in_band, "in_band")
  matrix$expected_lower <- as.numeric(matrix$expected_lower)
  matrix$expected_upper <- as.numeric(matrix$expected_upper)
  matrix$coverage <- as.numeric(matrix$coverage)
  matrix$n_replications_expected <- as.integer(matrix$n_replications_expected)
  matrix$n_replications_fit <- as.integer(matrix$n_replications_fit)
  matrix$n_replications_failed <- as.integer(matrix$n_replications_failed)
  matrix$n_sites <- as.integer(matrix$n_sites)

  expected_ids <- tier3_target_id(tier3_release_k_values())
  if (!setequal(matrix$target_id, expected_ids) || nrow(matrix) != length(expected_ids)) {
    stop("Tier 3 full-live matrix must contain exactly the release K targets.", call. = FALSE)
  }
  if (anyDuplicated(matrix$target_id) || anyDuplicated(matrix$K)) {
    stop("Tier 3 full-live matrix contains duplicate targets or K values.", call. = FALSE)
  }
  matrix <- matrix[match(expected_ids, matrix$target_id), , drop = FALSE]

  for (i in seq_len(nrow(matrix))) {
    row <- matrix[i, , drop = FALSE]
    target <- .bef_target(row$target_id, targets = targets, statuses = "active")
    if (!identical(row$K, as.integer(sub(".*_K", "", row$target_id)))) {
      stop("Tier 3 full-live matrix K does not match target_id.", call. = FALSE)
    }
    if (!identical(row$status, target$status)) {
      stop("Tier 3 full-live matrix status does not match active ledger.", call. = FALSE)
    }
    if (!identical(row$release_blocking, target$release_blocking)) {
      stop("Tier 3 full-live matrix release_blocking does not match ledger.", call. = FALSE)
    }
    if (!identical(row$fixture_path, target$fixture_path)) {
      stop("Tier 3 full-live matrix fixture_path does not match ledger.", call. = FALSE)
    }
    if (!isTRUE(all.equal(row$expected_lower, target$expected_value_lower))) {
      stop("Tier 3 full-live matrix lower bound does not match ledger.", call. = FALSE)
    }
    if (!isTRUE(all.equal(row$expected_upper, target$expected_value_upper))) {
      stop("Tier 3 full-live matrix upper bound does not match ledger.", call. = FALSE)
    }
    if (!isTRUE(row$in_band)) {
      stop("Tier 3 full-live matrix includes out-of-band coverage.", call. = FALSE)
    }
    if (!identical(row$n_replications_expected, 20L) ||
      !identical(row$n_replications_fit, 20L) ||
      !identical(row$n_replications_failed, 0L)) {
      stop("Tier 3 full-live matrix replication accounting is invalid.", call. = FALSE)
    }
    expected_sites <- 20L * as.integer(row$K)
    if (!identical(row$n_sites, expected_sites)) {
      stop("Tier 3 full-live matrix site denominator is invalid.", call. = FALSE)
    }
    if (!is.finite(row$coverage) ||
      row$coverage < target$expected_value_lower ||
      row$coverage > target$expected_value_upper) {
      stop("Tier 3 full-live matrix coverage is outside the active ledger bounds.", call. = FALSE)
    }
    if (!nzchar(row$evidence_path) || !dir.exists(row$evidence_path)) {
      stop("Tier 3 full-live matrix evidence_path does not exist.", call. = FALSE)
    }
  }

  row.names(matrix) <- NULL
  matrix
}

tier3_full_live_matrix_row <- function(matrix, K) {
  row <- matrix[matrix$K == as.integer(K), , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Tier 3 full-live matrix row lookup failed.", call. = FALSE)
  }
  row
}

tier3_target_active <- function(K, targets = .bef_load_targets()) {
  target <- .bef_target(tier3_target_id(K), targets = targets)
  identical(target$status, "active")
}

tier3_fit_config <- function(full = live_cmdstan_full_run_requested()) {
  if (isTRUE(full)) {
    return(list(chains = 4L, iter_warmup = 1000L, iter_sampling = 3000L))
  }
  list(chains = 1L, iter_warmup = 150L, iter_sampling = 20L)
}

with_live_cmdstan_cache_root <- function(code) {
  old <- Sys.getenv("BAYESEFRON_CACHE_ROOT", unset = NA_character_)
  root <- tempfile("bayesefron-live-cmdstan-cache-")
  reset_cache <- getFromNamespace(".bef_reset_session_cache", "bayesEfron")

  reset_cache()
  Sys.setenv(BAYESEFRON_CACHE_ROOT = root)
  on.exit({
    reset_cache()
    unlink(root, recursive = TRUE, force = TRUE)
    if (is.na(old)) {
      Sys.unsetenv("BAYESEFRON_CACHE_ROOT")
    } else {
      Sys.setenv(BAYESEFRON_CACHE_ROOT = old)
    }
  }, add = TRUE)

  force(code)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
