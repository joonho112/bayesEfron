tier3_schema_fixture <- function(K) {
  controls <- tier3_fixture_controls(K)
  K <- as.integer(K)
  replication <- tier3_schema_replication(1L, K, controls)

  list(
    metadata = tier3_schema_metadata(K, controls, replication, version = "v1"),
    K = replication$K,
    replication = replication$replication,
    dataset_id = replication$dataset_id,
    seed = replication$seed,
    theta_hat = replication$theta_hat,
    sigma = replication$sigma,
    theta_true = replication$theta_true,
    source_row_id = replication$source_row_id,
    tower = replication$tower,
    L = replication$L,
    M = replication$M,
    grid_method = replication$grid_method,
    expansion = replication$expansion,
    bound_expansion = replication$bound_expansion
  )
}

tier3_schema_fixture_v2 <- function(K) {
  controls <- tier3_fixture_controls(K)
  K <- as.integer(K)
  replications <- lapply(seq_len(20L), tier3_schema_replication, K = K, controls = controls)
  primary <- replications[[1L]]

  list(
    metadata = tier3_schema_metadata(K, controls, primary, version = "v2"),
    K = primary$K,
    replication = primary$replication,
    dataset_id = primary$dataset_id,
    seed = primary$seed,
    theta_hat = primary$theta_hat,
    sigma = primary$sigma,
    theta_true = primary$theta_true,
    source_row_id = primary$source_row_id,
    tower = primary$tower,
    L = primary$L,
    M = primary$M,
    grid_method = primary$grid_method,
    expansion = primary$expansion,
    bound_expansion = primary$bound_expansion,
    replications = replications
  )
}

tier3_schema_replication <- function(replication, K, controls) {
  replication <- as.integer(replication)
  K <- as.integer(K)
  left <- controls$left_count
  right <- controls$right_count

  theta_true <- c(seq(-2, -1, length.out = left), seq(1, 2, length.out = right))
  theta_hat <- theta_true + seq(-0.2, 0.2, length.out = K) + (replication - 1L) * 0.001
  sigma <- rep(0.5, K)
  source_row_id <- seq_len(K)
  tower <- c(rep("left_theta_true_lt_0", left), rep("right_theta_true_ge_0", right))
  observed_range <- range(theta_hat)
  grid_expansion <- diff(observed_range) * 0.5

  list(
    replication = replication,
    dataset_id = controls$dataset_id + replication - 1L,
    seed = controls$seed + replication - 1L,
    K = K,
    theta_hat = as.numeric(theta_hat),
    sigma = as.numeric(sigma),
    theta_true = as.numeric(theta_true),
    source_row_id = as.integer(source_row_id),
    tower = tower,
    L = controls$L,
    M = 6L,
    grid_method = "paper_realdata",
    expansion = 0.5,
    bound_expansion = 0.5,
    grid = list(
      observed_range = as.numeric(observed_range),
      grid_bounds = as.numeric(c(
        min(theta_hat) - grid_expansion,
        max(theta_hat) + grid_expansion
      )),
      grid_expansion = as.numeric(grid_expansion)
    )
  )
}

tier3_schema_metadata <- function(K, controls, primary, version = c("v1", "v2")) {
  version <- match.arg(version)
  metadata <- list(
    fixture_format_version = version,
    generator = "tools/generate-tier3-lee-sui-fixtures.R",
    generated_on = "2026-05-13",
    family = "Lee-Sui public small-K twin-towers",
    source = list(
      scenario = list(I = 0.7, R = 9),
      canonical_input = if (identical(as.integer(K), 1500L)) {
        "part01_full_panel_simulation_rule"
      } else {
        "archived_appendix_small_k_result"
      },
      archived_control_input = "archived_appendix_small_k_result",
      artifacts = tier3_schema_artifacts(),
      source_row_ids = primary$source_row_id,
      first_source_row_ids = utils::head(primary$source_row_id, 8L)
    ),
    K_scenarios = tier3_release_k_values(),
    n_replications = 20L,
    replication = 1L,
    dataset_id = controls$dataset_id,
    seed_base = 2025L,
    seed_range = c(2026L, 2125L),
    sampling_strategy = if (identical(as.integer(K), 1500L)) {
      "part01_seeded_full_panel_resimulation"
    } else {
      "theta_true_split_at_0_proportional_towers_unsorted"
    },
    RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
    grid = list(
      L = controls$L,
      M = 6L,
      grid_method = "paper_realdata",
      expansion = 0.5,
      bound_expansion = 0.5,
      observed_range = primary$grid$observed_range,
      grid_bounds = primary$grid$grid_bounds,
      grid_expansion = primary$grid$grid_expansion
    )
  )

  if (identical(version, "v2")) {
    metadata$K <- as.integer(K)
    metadata$primary_replication <- 1L
    metadata$replications <- seq_len(20L)
    metadata$dataset_id_range <- c(controls$dataset_id, controls$dataset_id_end)
    metadata$replication_seed_range <- c(controls$seed, controls$seed_end)
    metadata$grid$per_replication_observed_range <- TRUE
    metadata$coverage <- list(
      nominal_level = 0.9,
      aggregation = "mean_of_replication_site_coverages",
      n_replications = 20L,
      denominator = 20L * as.integer(K),
      failure_policy = "do_not_silently_drop_failed_replications",
      replication_identity_policy = "unique_input_payload_hashes",
      sampler_seed_policy = "tier3_sampler_seed"
    )
  }

  metadata
}

tier3_schema_artifacts <- function() {
  digest <- paste(rep("a", 64L), collapse = "")
  path <- "/tmp/source-artifact"
  artifact <- function(role, relevant_lines = "all") {
    list(
      role = role,
      path = path,
      sha256 = digest,
      relevant_lines = relevant_lines
    )
  }

  list(
    data_generation_script = artifact("public_data_generation_script", "65-156,327-331"),
    model_script = artifact("public_model_performance_script", "73-117,188-214"),
    public_small_k_script = artifact("public_small_k_script", "69-152,194-228,302-338,400-489"),
    original_small_k_script = artifact("original_small_k_development_script", "33-100,122-155,215-240,704-705"),
    stan_model = artifact("public_stan_model"),
    simulation_rds = artifact("public_simulation_rds", NA_character_),
    archived_small_k_result = artifact("archived_appendix_small_k_result", NA_character_)
  )
}

test_that("Tier 3 fixture controls encode the Lee-Sui K-series contract", {
  controls <- do.call(rbind, lapply(tier3_release_k_values(), tier3_fixture_controls))

  expect_equal(controls$K, tier3_release_k_values())
  expect_equal(controls$L, c(51L, 71L, 81L, 101L, 101L))
  expect_equal(controls$dataset_id, c(1L, 21L, 41L, 61L, 81L))
  expect_equal(controls$seed, c(2026L, 2046L, 2066L, 2086L, 2106L))
  expect_equal(controls$dataset_id_end, c(20L, 40L, 60L, 80L, 100L))
  expect_equal(controls$seed_end, c(2045L, 2065L, 2085L, 2105L, 2125L))
  expect_equal(controls$left_count + controls$right_count, controls$K)
})

test_that("Tier 3 fixture replication controls encode all 20 Lee-Sui replications", {
  for (K in tier3_release_k_values()) {
    controls <- tier3_fixture_controls(K)
    rep_controls <- tier3_replication_controls(K)

    expect_equal(rep_controls$replication, seq_len(20L))
    expect_equal(rep_controls$dataset_id, seq.int(controls$dataset_id, controls$dataset_id_end))
    expect_equal(rep_controls$seed, seq.int(controls$seed, controls$seed_end))
    expect_equal(rep_controls$seed, 2025L + rep_controls$dataset_id)
  }
})

test_that("Tier 3 fixture validator accepts v1 structure-only fixtures", {
  for (K in tier3_release_k_values()) {
    fixture <- tier3_schema_fixture(K)
    tier3_validate_fixture(fixture, K)
  }
})

test_that("Tier 3 fixture validator accepts v2 all-replication fixtures", {
  for (K in tier3_release_k_values()) {
    fixture <- tier3_schema_fixture_v2(K)
    tier3_validate_fixture(fixture, K)
  }
})

test_that("Tier 3 Lee-Sui fixtures on disk conform to schema", {
  fixture_files <- list.files(
    .bef_testthat_file("_fixtures"),
    pattern = "^lee_sui_K[0-9]+[.]rds$"
  )
  K_values <- sort(as.integer(sub("^lee_sui_K([0-9]+)[.]rds$", "\\1", fixture_files)))

  expect_setequal(K_values, tier3_release_k_values())

  for (K in K_values) {
    fixture <- tier3_load_fixture(K)
    expect_equal(fixture$metadata$fixture_format_version, "v2")
    expect_length(fixture$replications, 20L)
    tier3_validate_fixture(fixture, K)
  }
})

test_that("Tier 3 fixture validator rejects malformed v1 structure", {
  fixture <- tier3_schema_fixture(50L)
  fixture$metadata$source$canonical_input <- "public_part06_replay"

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    "archived_appendix_small_k_result"
  )

  fixture <- tier3_schema_fixture(50L)
  fixture$source_row_id[1L] <- fixture$source_row_id[2L]
  fixture$metadata$source$source_row_ids <- fixture$source_row_id
  fixture$metadata$source$first_source_row_ids <- utils::head(fixture$source_row_id, 8L)

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    "anyDuplicated|0"
  )

  fixture <- tier3_schema_fixture(50L)
  fixture$L <- 101L

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    "controls"
  )
})

test_that("Tier 3 fixture validator rejects malformed v2 structure", {
  fixture <- tier3_schema_fixture_v2(50L)
  fixture$replications <- NULL

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  fixture$replications <- fixture$replications[-20L]

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  fixture$theta_hat[1L] <- fixture$theta_hat[1L] + 1

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  fixture$replications[[2L]]$source_row_id[1L] <- fixture$replications[[2L]]$source_row_id[2L]

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  fixture$replications[[2L]]$seed <- fixture$replications[[2L]]$seed + 1L

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  payload_fields <- c("theta_hat", "sigma", "theta_true", "source_row_id", "tower", "grid")
  fixture$replications[[2L]][payload_fields] <- fixture$replications[[1L]][payload_fields]

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(50L)
  fixture$metadata$coverage$denominator <- 50L

  expect_error(
    tier3_validate_fixture(fixture, 50L),
    class = "expectation_failure"
  )

  fixture <- tier3_schema_fixture_v2(1500L)
  fixture$replications[[2L]]$source_row_id <- rev(fixture$replications[[2L]]$source_row_id)

  expect_error(
    tier3_validate_fixture(fixture, 1500L),
    class = "expectation_failure"
  )
})
