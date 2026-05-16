bef_test_source_root <- function(required = c("DESCRIPTION", "NAMESPACE")) {
  candidate_starts <- c(
    Sys.getenv("BAYESEFRON_SOURCE_ROOT", unset = NA_character_),
    getwd(),
    testthat::test_path()
  )
  candidate_starts <- candidate_starts[!is.na(candidate_starts)]
  parents <- unique(unlist(lapply(candidate_starts, bef_test_parent_chain)))
  rcheck_dirs <- parents[grepl("[.]Rcheck$", basename(parents))]
  candidates <- unique(c(
    parents,
    file.path(rcheck_dirs, "00_pkg_src", "bayesEfron")
  ))
  for (path in candidates) {
    if (all(file.exists(file.path(path, required)))) {
      return(normalizePath(path, mustWork = TRUE))
    }
  }
  stop("Cannot locate bayesEfron source root for source-file test.", call. = FALSE)
}

bef_test_parent_chain <- function(path, depth = 8L) {
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
