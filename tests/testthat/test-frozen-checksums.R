find_locked_source_root <- function() {
  candidate_starts <- c(
    Sys.getenv("BAYESEFRON_SOURCE_ROOT", unset = NA_character_),
    getwd(),
    testthat::test_path()
  )
  candidate_starts <- candidate_starts[!is.na(candidate_starts)]

  parents <- unique(unlist(lapply(candidate_starts, parent_chain)))
  rcheck_dirs <- parents[grepl("[.]Rcheck$", basename(parents))]
  candidates <- unique(c(
    parents,
    file.path(rcheck_dirs, "00_pkg_src", "bayesEfron")
  ))
  for (path in candidates) {
    if (is_locked_source_root(path)) {
      return(normalizePath(path, mustWork = TRUE))
    }
  }

  stop(
    "Cannot locate bayesEfron source root for locked-core checksum test.",
    call. = FALSE
  )
}

parent_chain <- function(path, depth = 8L) {
  out <- character()
  current <- normalizePath(path, mustWork = FALSE)
  for (i in seq_len(depth)) {
    out <- c(out, current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  out
}

is_locked_source_root <- function(path) {
  file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, "inst", "locked-core-checksums.txt")) &&
    file.exists(file.path(path, "inst", "stan", "efron_re.stan")) &&
    file.exists(file.path(path, "R", "grid.R")) &&
    file.exists(file.path(path, "R", "data-prep.R"))
}

test_that("locked core files match committed SHA-256 digests", {
  source_root <- find_locked_source_root()
  checksums_file <- file.path(source_root, "inst", "locked-core-checksums.txt")
  expect_true(file.exists(checksums_file))

  records <- readLines(checksums_file, warn = FALSE)
  expect_length(records, 3L)

  parts <- strsplit(records, "  ", fixed = TRUE)
  expect_true(all(lengths(parts) == 2L))

  expected <- stats::setNames(
    vapply(parts, `[[`, character(1), 1L),
    vapply(parts, `[[`, character(1), 2L)
  )
  expected_paths <- c(
    "inst/stan/efron_re.stan",
    "R/grid.R",
    "R/data-prep.R"
  )

  expect_setequal(names(expected), expected_paths)
  expect_equal(anyDuplicated(names(expected)), 0L)
  expect_true(all(grepl("^[[:xdigit:]]{64}$", expected)))

  actual <- vapply(
    names(expected),
    function(path) {
      source_path <- file.path(source_root, path)
      expect_true(file.exists(source_path), info = path)
      digest::digest(
        source_path,
        algo = "sha256",
        file = TRUE
      )
    },
    character(1)
  )

  expect_equal(actual[names(expected)], expected)
})
