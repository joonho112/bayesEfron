locked_tier3_manifest_path <- function() {
  .bef_source_file(
    "inst",
    "locked-tier3-fixtures-checksums.txt",
    required = c("DESCRIPTION", "inst/locked-tier3-fixtures-checksums.txt")
  )
}

read_locked_tier3_manifest <- function(path = locked_tier3_manifest_path()) {
  utils::read.delim(
    path,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

test_that("locked Tier 3 Lee-Sui fixtures match committed SHA-256 digests", {
  manifest <- read_locked_tier3_manifest()
  K_values <- tier3_release_k_values()
  expected_ids <- sprintf("tier3_fixture_lee_sui_K%d_sha256", K_values)
  expected_paths <- sprintf("tests/testthat/_fixtures/lee_sui_K%d.rds", K_values)

  expect_named(manifest, c("target_id", "path", "sha256", "provenance"))
  expect_length(manifest$target_id, length(K_values))
  expect_equal(anyDuplicated(manifest$target_id), 0L)
  expect_equal(anyDuplicated(manifest$path), 0L)
  expect_setequal(manifest$target_id, expected_ids)
  expect_setequal(manifest$path, expected_paths)
  expect_true(all(grepl("^[[:xdigit:]]{64}$", manifest$sha256)))
  expect_true(all(nzchar(manifest$provenance)))
  expect_false(any(grepl("^/|(^|/)[.][.](/|$)", manifest$path)))

  targets <- .bef_load_targets()
  release_rows <- targets[match(tier3_target_id(K_values), targets$target_id), , drop = FALSE]
  expect_false(anyNA(release_rows$target_id))
  expect_setequal(manifest$path, file.path("tests", "testthat", release_rows$fixture_path))

  actual <- vapply(
    manifest$path,
    function(path) {
      source_path <- .bef_source_file(path, required = c("DESCRIPTION", path))
      expect_true(file.exists(source_path), info = path)
      digest::digest(source_path, algo = "sha256", file = TRUE)
    },
    character(1)
  )
  names(actual) <- manifest$target_id

  expected <- stats::setNames(manifest$sha256, manifest$target_id)
  expect_equal(actual[names(expected)], expected)
})
