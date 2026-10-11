# Builds inst/examples/raudenbush_ranges.rds: compare_efron_grids() for the
# 19 teacher expectancy studies on three ranges, which the examples of
# ?compare_efron_grids and ?plot.bef_grid_comparison read. Needs CmdStan and
# the loo package. Run from the package root:
#
#   Rscript data-raw/range-comparison.R
#
# The three ranges widen the range of the estimates on each side by 0.25,
# 0.5 and 1 times its width (rule "paper_realdata"), with 101 grid points
# and six spline functions. The comparison keeps the summaries of studies 4
# and 8, the studies with the largest and the smallest estimate. Each fit
# runs four chains of 1,000 warmup and 3,000 sampling iterations.

suppressMessages(devtools::load_all(".", quiet = TRUE))

raudenbush_ranges <- suppressMessages(
  compare_efron_grids(raudenbush1985, L = 101, M = 6, expansion = c(0.25, 0.5, 1),
                      sites = c(4, 8), seed = 1985),
  classes = "bef_progress"
)
stopifnot(
  nrow(raudenbush_ranges) == 3L, !anyNA(raudenbush_ranges$elpd_loo),
  identical(which.max(raudenbush1985$yi), 4L), identical(which.min(raudenbush1985$yi), 8L)
)
# How long the sampler took is not part of what the file shows.
raudenbush_ranges$sampler_seconds <- NA_real_
destination <- file.path("inst", "examples", "raudenbush_ranges.rds")
saveRDS(raudenbush_ranges, destination, compress = "xz", version = 3)
print(raudenbush_ranges)
message(sprintf("%s: %d bytes", destination, file.size(destination)))
