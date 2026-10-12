# ==============================================================================
# ratio-mcse.R: paired Monte Carlo standard errors for the length and RMSE
# ratios in Figure 2 (the plotting file is named fig4_evidence).
#
# Run from the paper folder root:
#   Rscript scripts/ratio-mcse.R /path/to/bayesEfron-checkout
# The checkout argument is optional when paper/ is inside the checkout. The
# historical sibling checkout names are also recognized for the author folder.
# Add --check-inputs to validate the checkout without producing output.
# Reads saved replication results only; it does not refit a model. The two-group
# panels are shipped with bayesEfron 0.4.0; the 19-study and equal-SE rows are in
# data-raw/benchmark/results-v040/ in the package source checkout.
#
# For each method X and the bayesEfron reference Y, the plotted ratio is
# r = mean(X) / mean(Y). Its paired delta-method Monte Carlo SE is
# r * sd(X / mean(X) - Y / mean(Y)) / sqrt(B), where B is the number of
# independently generated data sets. The same replication identifiers must
# occur once in each method. The reference divided by itself has zero SE.
#
# Writes data/ratio-mcse.csv, which lets figures.R draw the intervals without
# access to the package source checkout. Every mean is checked against the
# corresponding stored summary in data/evidence.rds before the CSV is written.
# ==============================================================================

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "scripts", "ratio-mcse.R")))
.libPaths(c(file.path(root, "r-library"), .libPaths()))
print(packageVersion("bayesEfron"))
stopifnot(packageVersion("bayesEfron") == "0.4.0")

args <- commandArgs(trailingOnly = TRUE)
check_inputs <- "--check-inputs" %in% args
args <- args[args != "--check-inputs"]
if (length(args) > 1L) {
  stop("Usage: Rscript scripts/ratio-mcse.R [package-checkout] [--check-inputs]",
       call. = FALSE)
}
candidates <- if (length(args)) args else c(
  file.path(root, ".."),
  file.path(root, "..", "bayesEfron-R-package-public"),
  file.path(root, "..", "bayesEfron-R-package")
)
needed <- c("few-studies-rows.csv", "equal-se-M6.csv")
has_results <- vapply(candidates, function(path) {
  all(file.exists(file.path(path, "data-raw", "benchmark", "results-v040", needed)))
}, logical(1))
if (!any(has_results)) {
  missing <- vapply(candidates, function(path) {
    absent <- needed[!file.exists(file.path(path, "data-raw", "benchmark",
                                           "results-v040", needed))]
    paste0("  ", path, ": missing ", paste(absent, collapse = ", "))
  }, character(1))
  stop("Cannot locate results-v040. Supply the bayesEfron source checkout: ",
       "Rscript scripts/ratio-mcse.R /path/to/bayesEfron-checkout\n",
       paste(missing, collapse = "\n"), call. = FALSE)
}
checkout <- normalizePath(candidates[which(has_results)[1L]])
bench <- file.path(checkout, "data-raw", "benchmark", "results-v040")
message("benchmark folder: ", bench)
if (check_inputs) {
  message("All benchmark input files found; no results were written.")
  quit(save = "no", status = 0L)
}

evidence <- readRDS(file.path(root, "data", "evidence.rds"))
panels <- readRDS(system.file("examples", "ebnm_comparison.rds",
                            package = "bayesEfron"))$panels$by_panel
few <- read.csv(file.path(bench, "few-studies-rows.csv"), stringsAsFactors = FALSE)
equal <- read.csv(file.path(bench, "equal-se-M6.csv"), stringsAsFactors = FALSE)

designs <- data.frame(
  design_id = c("two_groups_K100", "two_groups_K1500", "few_normal_014",
                "few_normal_031", "few_fitted_g", "equal_se_K1500"),
  design = c("Two groups, K = 100\n(100 data sets)",
             "Two groups, K = 1,500\n(20 data sets)",
             "19 studies, normal SD 0.14\n(200 data sets)",
             "19 studies, normal SD 0.31\n(200 data sets)",
             "19 studies, fitted g\n(200 data sets)",
             "Two groups, equal SE, K = 1,500\n(20 data sets)"),
  replications = c(100L, 20L, 200L, 200L, 200L, 20L),
  stringsAsFactors = FALSE
)
panel_methods <- c("bayesEfron" = "bayesEfron",
                   "g of bayesEfron, plugged in" = "bayesEfron g, plugged in",
                   "Normal" = "Normal (plug-in)", "NPMLE" = "NPMLE")
few_methods <- c("bayesEfron" = "bayesEfron", "normal" = "Normal (plug-in)",
                 "known g" = "Known g")
equal_methods <- c("bayesEfron" = "bayesEfron",
                   "g of bayesEfron, plugged in" = "bayesEfron g, plugged in",
                   "Efron's estimator" = "Efron's estimator")

assert_near <- function(actual, expected, what) {
  if (length(actual) != 1L || length(expected) != 1L ||
      !is.finite(actual) || !is.finite(expected) || abs(actual - expected) > 1e-12) {
    stop("Summary mismatch: ", what, call. = FALSE)
  }
}
output <- diagnostics <- list()
add_design <- function(index, raw, summary, id_column, length_column, methods) {
  reference <- raw[raw$method == "bayesEfron", , drop = FALSE]
  reference <- reference[order(reference[[id_column]]), , drop = FALSE]
  stopifnot(nrow(reference) == designs$replications[index],
            !anyNA(reference[[id_column]]), !anyDuplicated(reference[[id_column]]))
  reference_summary <- summary[summary$method == "bayesEfron", , drop = FALSE]
  stopifnot(nrow(reference_summary) == 1L)
  for (method in names(methods)) {
    current <- raw[raw$method == method, , drop = FALSE]
    current <- current[order(current[[id_column]]), , drop = FALSE]
    stopifnot(!anyDuplicated(current[[id_column]]),
              identical(current[[id_column]], reference[[id_column]]))
    current_summary <- summary[summary$method == method, , drop = FALSE]
    stopifnot(nrow(current_summary) == 1L)
    row <- data.frame(design_id = designs$design_id[index],
                      design = designs$design[index], method = unname(methods[method]),
                      replications = nrow(current), stringsAsFactors = FALSE)
    for (metric in c("length", "rmse")) {
      column <- if (metric == "length") length_column else "rmse"
      x <- current[[column]]
      y <- reference[[column]]
      stopifnot(all(is.finite(x)), all(is.finite(y)), mean(x) > 0, mean(y) > 0)
      assert_near(mean(x), current_summary[[column]],
                  paste(designs$design_id[index], method, metric))
      assert_near(mean(y), reference_summary[[column]],
                  paste(designs$design_id[index], "bayesEfron", metric))
      ratio <- mean(x) / mean(y)
      assert_near(ratio, current_summary[[column]] / reference_summary[[column]],
                  paste(designs$design_id[index], method, metric, "ratio"))
      se <- if (method == "bayesEfron") 0 else
        ratio * stats::sd(x / mean(x) - y / mean(y)) / sqrt(length(x))
      row[[paste0(metric, "_ratio")]] <- ratio
      row[[paste0(metric, "_ratio_se")]] <- se
      diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
        design_id = designs$design_id[index], method = unname(methods[method]),
        metric = metric, covariance = stats::cov(x, y), ratio = ratio,
        unpaired_se = ratio * sqrt((stats::sd(x) / mean(x))^2 +
                                  (stats::sd(y) / mean(y))^2) / sqrt(length(x)),
        paired_se = se, stringsAsFactors = FALSE
      )
    }
    output[[length(output) + 1L]] <<- row
  }
}

for (K in c(100L, 1500L)) {
  add_design(if (K == 100L) 1L else 2L, panels[panels$K == K, ],
             evidence$panels_methods[evidence$panels_methods$K == K, ],
             "panel", "mean_length", panel_methods)
}
for (i in seq_along(scenarios <- c("normal, SD 0.14", "normal, SD 0.31", "fitted g"))) {
  add_design(i + 2L, few[few$scenario == scenarios[i], ],
             evidence$few[evidence$few$scenario == scenarios[i], ],
             "replication", "length", few_methods)
}
add_design(6L, equal[equal$M == 6, ], evidence$equal[evidence$equal$M == 6, ],
           "panel", "mean_length", equal_methods)

output <- do.call(rbind, output)
diagnostics <- do.call(rbind, diagnostics)
stopifnot(nrow(output) == 20L,
          !anyDuplicated(output[c("design_id", "method")]),
          all(output$length_ratio_se >= 0), all(output$rmse_ratio_se >= 0),
          all(output$length_ratio_se[output$method == "bayesEfron"] == 0),
          all(output$rmse_ratio_se[output$method == "bayesEfron"] == 0))
counterexample <- diagnostics[diagnostics$design_id == "two_groups_K1500" &
                               diagnostics$method == "Normal (plug-in)" &
                               diagnostics$metric == "length", ]
stopifnot(nrow(counterexample) == 1L, counterexample$covariance < 0,
          counterexample$paired_se > counterexample$unpaired_se)
write.csv(output, file.path(root, "data", "ratio-mcse.csv"), row.names = FALSE)
print(output, digits = 9, row.names = FALSE)
cat("\nNegative covariance example (paired SE exceeds the unpaired approximation):\n")
print(counterexample, digits = 10, row.names = FALSE)
message("All 40 point ratios match the existing evidence summaries within 1e-12.")
message("Wrote data/ratio-mcse.csv: 20 design-method rows, with paired SEs for both ratios.")
