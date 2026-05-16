locked_input_manifest_path <- function() {
  .bef_source_file(
    "inst",
    "locked-input-fixtures-checksums.txt",
    required = c("DESCRIPTION", "inst/locked-input-fixtures-checksums.txt")
  )
}

read_locked_input_manifest <- function(path = locked_input_manifest_path()) {
  utils::read.delim(
    path,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

test_that("locked input fixtures match committed SHA-256 digests", {
  manifest <- read_locked_input_manifest()

  expect_named(manifest, c("target_id", "path", "sha256", "provenance"))
  expect_length(manifest$target_id, 3L)
  expect_equal(anyDuplicated(manifest$target_id), 0L)
  expect_equal(anyDuplicated(manifest$path), 0L)

  expect_setequal(
    manifest$target_id,
    c(
      "tier0_input_fixture_paper_realdata_sha256",
      "tier0_input_fixture_paper_simulation_sha256",
      "tier0_input_fixture_paper_sensitivity_sha256"
    )
  )
  expect_setequal(
    manifest$path,
    c(
      "tests/testthat/_fixtures/grid/paper_realdata_inputs.rds",
      "tests/testthat/_fixtures/grid/paper_simulation_inputs.rds",
      "tests/testthat/_fixtures/grid/paper_sensitivity_inputs.rds"
    )
  )
  expect_true(all(grepl("^[[:xdigit:]]{64}$", manifest$sha256)))
  expect_true(all(nzchar(manifest$provenance)))
  expect_false(any(grepl("^/|(^|/)[.][.](/|$)", manifest$path)))

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
