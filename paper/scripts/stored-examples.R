# ==============================================================================
# stored-examples.R: numbers that the supplement quotes from the package's
# stored examples (inst/examples/lee_sui_examples.rds): the sensitivity of
# the count of firms with P(gap > 0) >= 0.95 to the grid, and the share of g
# above 5 percentage points in each stored estimate of g.
#
# Run from the paper folder root:  Rscript scripts/stored-examples.R
# Writes data/stored-examples-numbers.csv.
# ==============================================================================
root <- normalizePath(".")
.libPaths(c(file.path(root, "r-library"), .libPaths()))
path <- system.file("examples", "lee_sui_examples.rds", package = "bayesEfron")
stopifnot(nzchar(path))
ex <- readRDS(path)
sens <- ex$sensitivity[ex$sensitivity$case == "krw", ]
rownames(sens) <- sens$spec
reg <- list(
  krw_stored_positive_default = sens["baseline", "positive_intervals"],
  krw_stored_positive_narrower = sens["narrower", "positive_intervals"],
  krw_stored_positive_basis8 = sens["basis8", "positive_intervals"],
  krw_stored_sd_range_pp = 100 * diff(range(sens$sd_g))
)
for (spec in c("baseline", "narrower", "basis8")) {
  cur <- ex$curves[ex$curves$case == "krw" & ex$curves$spec == spec, ]
  stopifnot(abs(sum(cur$probability) - 1) < 1e-6)
  reg[[paste0("krw_stored_above5_", spec)]] <- sum(cur$probability[cur$grid > 0.05])
}
out <- data.frame(name = names(reg), value = unlist(reg), stringsAsFactors = FALSE)
write.csv(out, "data/stored-examples-numbers.csv", row.names = FALSE)
print(out, digits = 4)
