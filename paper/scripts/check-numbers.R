# ==============================================================================
# check-numbers.R: checks registered numerical quotations in main.tex and the
# supplement against the preserved registries in data/reference/.
#
# Run from the paper folder root with --mode=archival or --mode=reproduction.
# Archival mode requires manuscript sources and checks their recorded quotations
# against preserved references, including the recorded timings. Reproduction mode
# compares current and preserved registries without manuscript sources; timings
# and environment versions are reported separately from statistical differences.
# Neither mode overwrites the reference registries.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  matched <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(matched) > 1L) stop("Duplicate option: ", name)
  if (length(matched)) sub(paste0("^--", name, "="), "", matched) else default
}
if (any(!grepl("^--(mode|reference)=", args))) {
  stop("Usage: Rscript scripts/check-numbers.R --mode=archival|reproduction --reference=data/reference")
}
mode <- option("mode", "archival")
if (!mode %in% c("archival", "reproduction")) stop("Unknown mode: ", mode)
reference <- option("reference", "data/reference")
registries <- c("numbers.csv", "stored-examples-numbers.csv", "tail-coverage-numbers.csv")
read_registry <- function(path) {
  if (!file.exists(path)) stop("Missing registry: ", path)
  x <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
  if (!all(c("name", "value") %in% names(x))) stop("Registry needs name and value columns: ", path)
  if (anyNA(x$name) || any(!nzchar(x$name)) || anyDuplicated(x$name)) {
    stop("Missing or duplicate registry key: ", path)
  }
  if (anyNA(x$value) || any(!nzchar(x$value))) stop("Missing registry value: ", path)
  x
}
.libPaths(c(file.path(getwd(), "r-library"), .libPaths()))
print(packageVersion("bayesEfron"))
stopifnot(packageVersion("bayesEfron") == "0.4.0")
if (mode == "reproduction") {
  failures <- character()
  timing <- function(x) grepl("(^|_)(seconds|minutes)(_|$)", x)
  environment <- c("r_version", "cmdstan_version", "cmdstanr_version", "ebnm_version",
                   "metafor_version", "platform", "analysis_date")
  for (name in registries) {
    paths <- c(file.path(reference, name), file.path("data", name))
    old <- read_registry(paths[1])
    new <- read_registry(paths[2])
    if (!setequal(old$name, new$name)) stop("Registry keys differ: ", name)
    new <- new[match(old$name, new$name), ]
    for (i in seq_len(nrow(old))) {
      key <- old$name[i]; a <- old$value[i]; b <- new$value[i]
      if (timing(key) || key %in% environment) {
        cat(if (timing(key)) "TIMING" else "ENVIRONMENT", key, ": recorded=", a, "; rerun=", b, "\n", sep = " ")
        next
      }
      aa <- suppressWarnings(as.numeric(a)); bb <- suppressWarnings(as.numeric(b))
      ok <- if (is.finite(aa) && is.finite(bb)) abs(aa-bb) <= 1e-8 * max(1, abs(aa)) else identical(a,b)
      if (!ok) failures <- c(failures, sprintf("%s: %s recorded=%s rerun=%s",name,key,a,b))
    }
  }
  if (length(failures)) {
    cat("STATISTICAL DIFFERENCES (not timing failures):\n", paste(failures, collapse="\n"), "\n")
    cat("The 1e-8 tolerance allows rounding only. Different platforms or software may require a separate Monte Carlo error assessment.\n")
    quit(status=1L)
  }
  cat("REPRODUCTION-OK: all non-timing statistical registry values match the preserved reference (relative/absolute tolerance 1e-8).\n")
  quit(status=0L)
}
if (!file.exists("main.tex") || !dir.exists("osm/appendices")) {
  stop("Archival mode needs the manuscript sources supplied in the supplementary ZIP; use --mode=reproduction for a rerun comparison.")
}

cat("ARCHIVAL references:", reference, "(recorded quotations, including timings)\n")
numbers <- read_registry(file.path(reference, "numbers.csv"))
val <- setNames(numbers$value, numbers$name)
shown <- setNames(numbers$shown, numbers$name)
# the smaller registries of tail-coverage.R and stored-examples.R
for (fx in file.path(reference, c("tail-coverage-numbers.csv", "stored-examples-numbers.csv"))) {
  e <- read_registry(fx)
  val <- c(val[setdiff(names(val), e$name)], setNames(as.character(e$value), e$name))
}
main <- paste(readLines("main.tex", warn = FALSE), collapse = "\n")
osm <- paste(unlist(lapply(list.files("osm/appendices", full.names = TRUE), readLines, warn = FALSE)),
             collapse = "\n")
text <- paste(main, osm)
failures <- character()
squash <- function(x) gsub("\\s+", " ", x)
text <- squash(text)
expect <- function(what, pattern, where = text) {
  if (!grepl(squash(pattern), where, fixed = TRUE)) failures <<- c(failures, paste(what, "->", pattern))
}
num <- function(name) as.numeric(val[[name]])
f <- function(x, d) formatC(x, format = "f", digits = d)

# ---- the 97 firms --------------------------------------------------------------
expect("firms", sprintf("%s firms", val[["krw_K"]]))
expect("applications", format(num("krw_applications"), big.mark = ","))
expect("se ratio", sprintf("%s times the smallest", f(num("krw_se_ratio"), 1)))
expect("grid lower pp", sprintf("$-%s$ to $%s$ percentage points", f(-num("krw_grid_lower_pp"), 1),
                                 f(num("krw_grid_upper_pp"), 1)))
expect("mean g", sprintf("%s points (90\\%% posterior interval %s to", f(num("krw_mean_g_pp"), 2),
                         f(num("krw_mean_g_lower_pp"), 2)))
expect("mean g upper", sprintf("to\n%s)", f(num("krw_mean_g_upper_pp"), 2)))
expect("sd g", sprintf("is %s (%s to %s)", f(num("krw_sd_g_pp"), 2), f(num("krw_sd_g_lower_pp"), 2),
                       f(num("krw_sd_g_upper_pp"), 2)))
expect("sd estimates", sprintf("97 estimates, %s", f(num("krw_sd_estimates_pp"), 2)))
expect("sd posterior means", sprintf("posterior means, %s", f(num("krw_sd_posterior_means_pp"), 2)))
expect("reml mean", sprintf("same mean, %s", f(num("krw_reml_mean_pp"), 2)))
expect("reml sd", sprintf("standard deviation, %s", f(num("krw_reml_sd_pp"), 2)))
expect("above 5 normal", sprintf("probability %s above 5", f(num("krw_normal_above_5"), 3)))
expect("above 5 g", sprintf("puts\n%s (%s to %s)", f(num("krw_prob_above_5"), 3),
                            f(num("krw_prob_above_5_lower"), 3), f(num("krw_prob_above_5_upper"), 3)))
expect("count 0.95", sprintf("The %s firms that", val[["krw_n_prob_ge_95"]]))
expect("top estimate", sprintf("largest estimate, %s points", f(num("krw_top_estimate_pp"), 1)))
expect("top mean", sprintf("posterior mean of %s", f(num("krw_top_mean_pp"), 1)))
expect("share median", sprintf("The median is %s\\%% for the firms", f(num("krw_share_about_g_median_pct"), 1)))
expect("krw seconds, supplement", sprintf("took %s seconds of wall-clock time; the sampling itself took %s seconds and the pilot fit %s seconds",
                                          f(num("krw_seconds_total"), 1), f(num("krw_seconds_sampling"), 1), f(num("krw_seconds_pilot"), 1)))
expect("twin seconds, supplement", sprintf("took %s seconds of wall-clock time: %s seconds for the pilot fit on the wide range, %s seconds of sampling",
                                           f(num("twin_seconds_total"), 0), f(num("twin_seconds_pilot"), 0), f(num("twin_seconds_sampling"), 0)))
expect("share max", sprintf("reaches %s\\%% for the firm", f(num("krw_share_about_g_max_pct"), 0)))
expect("share max, impact", sprintf("reaching %s\\%%", f(num("krw_share_about_g_max_pct"), 0)))
expect("sensitivity counts", sprintf("%s with a narrower range and %s with 8 spline functions",
                                     val[["krw_stored_positive_narrower"]], val[["krw_stored_positive_basis8"]]))
expect("sensitivity counts, supplement", sprintf("gives %s and %s firms", val[["krw_stored_positive_narrower"]],
                                                 val[["krw_stored_positive_basis8"]]))
expect("stored above 5, default", sprintf("gives %s for the default fit", f(num("krw_stored_above5_baseline"), 3)))
expect("stored above 5, others", sprintf("is %s with the narrower range and %s with 8", f(num("krw_stored_above5_narrower"), 3),
                                         f(num("krw_stored_above5_basis8"), 3)))
expect("analysis minutes", "about two and a half minutes")
stopifnot(abs(num("analysis_minutes") - 2.5) < 0.5)
expect("length ratio plugin", sprintf("shorter by %s\\%% with 100 sites", f(num("extra_length_K100_pct"), 1)))
expect("krw seconds", "about five seconds")
stopifnot(abs(num("krw_seconds_total") - 5) < 1.5)

# ---- the simulated sites -------------------------------------------------------
expect("twin groups", sprintf("%s values between", val[["twin_low_n"]]))
expect("twin se range", sprintf("standard errors from %s to %s", f(num("twin_se_min"), 2), f(num("twin_se_max"), 2)))
expect("twin se above gap", sprintf("for %s\\%% of the sites", f(num("twin_se_above_gap_pct"), 0)))
expect("twin wide range", sprintf("from $-%s$ to $%s$ to $-%s$ to $%s$", f(-num("twin_wide_lower"), 1),
                                  f(num("twin_wide_upper"), 1), f(-num("twin_grid_lower"), 1), f(num("twin_grid_upper"), 1)))
expect("twin seconds", sprintf("took %s seconds, %s of them for the pilot", f(num("twin_seconds_total"), 0),
                               f(num("twin_seconds_pilot"), 0)))
expect("twin leaning site", sprintf("site %s, whose estimate", val[["twin_leaning_site"]]))
expect("twin leaning mass", sprintf("has %s\\%% of its posterior", f(num("twin_leaning_mass_pct"), 0)))
expect("twin between", sprintf("puts probability %s on the empty interval", f(num("twin_prob_between"), 2)))
expect("twin normal between", sprintf("model puts %s;", f(num("twin_normal_between"), 2)))
expect("twin moments", sprintf("%s and %s, match the\ntrue %s and %s", f(num("twin_mean_g"), 2), f(num("twin_sd_g"), 2),
                               f(num("twin_mean_truth"), 2), f(num("twin_sd_truth2"), 2)))
expect("twin site a mean", sprintf("%s for site %s", f(num("twin_site_a_mean"), 2), val[["twin_site_a"]]))
expect("twin inside", sprintf("contain %s\\%% of the true", f(num("twin_inside_90_pct"), 1)))
expect("twin rmse", sprintf("is %s, against %s for the estimates and %s", f(num("twin_rmse_posterior"), 3),
                            f(num("twin_rmse_estimates"), 3), f(num("twin_rmse_normal"), 3)))
expect("twin npmle", sprintf("intervals contain %s\\%%", f(num("twin_eb_inclusion_npmle_pct"), 1)))

# ---- the 19 studies ------------------------------------------------------------
expect("raud sd g range", sprintf("between %s and %s", f(num("raud_sd_g_min"), 2), f(num("raud_sd_g_max"), 2)))
expect("raud reml sd", sprintf("REML estimate of %s", f(num("raud_reml_sd"), 2)))
expect("raud site 4", sprintf("from %s on the\nnarrowest range to %s and %s", f(num("raud_site4_narrow"), 2),
                              f(num("raud_site4_default"), 2), f(num("raud_site4_wide"), 2)))
expect("raud seconds", "take seven seconds in all")
stopifnot(abs(num("raud_seconds_total") - 7) < 1.5)

# ---- the evidence -------------------------------------------------------------
expect("cov K100", sprintf("contained %s\\%% (Monte Carlo\nstandard error %s)", f(100 * num("cov_K100"), 1),
                           f(100 * num("cov_K100_se"), 1)))
expect("cov K1500", sprintf("and %s\\%% (%s) of the true", f(100 * num("cov_K1500"), 1), f(100 * num("cov_K1500_se"), 1)))
expect("cov few", sprintf("%s\\%% to %s\\%%, whereas", f(100 * num("cov_few_min"), 1), f(100 * num("cov_few_max"), 1)))
expect("cov few normal", sprintf("contained %s\\%% to %s\\%%", f(100 * num("cov_few_normal_min"), 0),
                                 f(100 * num("cov_few_normal_max"), 0)))
expect("cov equal", sprintf("equal standard errors, %s\\%% (%s) against %s\\%%", f(100 * num("cov_equal_M6"), 1),
                            f(100 * num("cov_equal_M6_se"), 1), f(100 * num("cov_equal_efron_M6"), 1)))
expect("cov npmle", sprintf("contained\n%s\\%% and %s\\%% of the true", f(100 * num("cov_npmle_K100"), 0),
                            f(100 * num("cov_npmle_K1500"), 0)))
ev <- readRDS("data/evidence.rds")$few
fb <- ev[ev$method == "bayesEfron", ]
factor_truth <- fb$sd_mean / fb$true_sd
stopifnot(abs(num("few_sd_factor_min") - min(factor_truth)) < 1e-12,
          abs(num("few_sd_factor_max") - max(factor_truth)) < 1e-12)
expect("sd factor (relative to the true SD)", sprintf("factor of %s to %s", f(num("few_sd_factor_min"), 1), f(num("few_sd_factor_max"), 1)))
expect("panels valid", sprintf("All %s\nfits with the defaults are valid", val[["panels_default_valid"]]))
expect("panels seconds", sprintf("median of %s seconds", f(num("panels_default_seconds"), 0)))
expect("v030 valid", sprintf("%s of the 20 fits are valid", val[["panels_v030_valid"]]))
expect("v030 seconds", sprintf("takes %s seconds", f(num("panels_v030_seconds"), 0)))

# ---- package facts ------------------------------------------------------------
expect("version", "0.4.0")
expect("exports", "Ten plot functions return")

# ---- subsets of the sites (tail-coverage.R) ------------------------------------
pct <- function(name) f(100 * num(name), 0)
expect("tail extreme", sprintf("it is %s\\%% and %s\\%% (100 and 1,500 sites) outside the sample 2.5th and 97.5th percentiles",
                               pct("tail_K100_be_extreme"), pct("tail_K1500_be_extreme")))
expect("tail middle", sprintf("%s\\%% and %s\\%% between the groups", pct("tail_K100_be_middle"),
                              pct("tail_K1500_be_middle")))
expect("tail outer", sprintf("and %s\\%% and %s\\%% for true effects", pct("tail_K100_be_outer"), pct("tail_K1500_be_outer")))
stopifnot(pct("tail_K100_true_outer") == pct("tail_K1500_true_outer"))
expect("tail outer, true g", sprintf("The true $g$ gives %s\\%% at those edges in both designs", pct("tail_K100_true_outer")))
stopifnot(pct("tail_K100_normal_outer") == pct("tail_K1500_normal_outer"))
expect("tail outer, normal", sprintf("The normal model gives %s\\%% there", pct("tail_K100_normal_outer")))
expect("tail plug-in gap", sprintf("at most %s points", f(num("tail_K100_plug_max_gap_pp"), 1)))
stopifnot(num("tail_K1500_plug_max_gap_pp") < 0.1)
expect("tail scope", sprintf("contained %s\\%% to %s\\%% of the effects at the edges", pct("tail_K1500_be_outer"),
                             pct("tail_K100_be_outer")))
expect("tail supplement, extreme", sprintf("contain %s\\%% and %s\\%% (100 and 1,500 sites) of the true effects of estimates outside those sample percentiles",
                                           pct("tail_K100_be_extreme"), pct("tail_K1500_be_extreme")))
expect("tail supplement, middle true g", sprintf("where the true $g$ gives %s\\%% and %s\\%%;", pct("tail_K100_true_middle"),
                                                 pct("tail_K1500_true_middle")))
expect("tail supplement, normal middle", sprintf("%s\\%% and %s\\%% between the groups", pct("tail_K100_normal_middle"),
                                                 pct("tail_K1500_normal_middle")))
stopifnot(num("n_plot_functions") == 10, num("n_exported") == 21)
expect("cmdstan", sprintf("CmdStan %s", val[["cmdstan_version"]]))
expect("R version", sprintf("R %s", val[["r_version"]]))

if (length(failures)) {
  cat("MISMATCHES:\n"); cat(paste(" -", failures), sep = "\n")
  quit(status = 1)
}
cat("CHECK-NUMBERS-OK:", nrow(numbers), "primary registry entries;", "registered numerical quotations checked against the archived values.\n")
