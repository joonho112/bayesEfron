test_that("cli soft dependency is contained in format-print.R", {
  root <- bef_test_source_root(required = c("DESCRIPTION", "R"))
  r_dir <- file.path(root, "R")
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  hits <- unlist(
    lapply(
      r_files,
      function(path) {
        lines <- readLines(path, warn = FALSE)
        if (any(grepl("cli::|requireNamespace\\(\"cli\"", lines))) {
          return(path)
        }
        character()
      }
    ),
    use.names = FALSE
  )
  hits <- if (length(hits) > 0L) file.path("R", basename(hits)) else character()

  expect_equal(unique(hits), "R/format-print.R")
})
