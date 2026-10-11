# Builds inst/verification/few-studies.csv: what a fit of this package and
# the normal random-effects model report for 19 studies when the truth is
# known. vignette("teacher-expectancy") reads the table. Needs CmdStan. Run
# from the package root:
#
#   Rscript data-raw/benchmark/few-studies.R [replications] [workers] [rows-file]
#   Rscript data-raw/benchmark/few-studies.R --analytic-only [rows-file] [table-file]
#
# The second form makes no fits. It computes the known-g rows of the third
# distribution again, in the rows file and in the table, and leaves the
# other rows and the date of the table as they are. Its default rows file
# is data-raw/benchmark/results-v040/few-studies-rows.csv and its default
# table inst/verification/few-studies.csv.
#
# Every data set has the 19 standard errors of raudenbush1985, and an
# estimate is its true effect plus normal noise with its standard error. The
# true effects are drawn from one of three distributions:
#
#   normal, SD 0.14   a normal distribution with the mean and the standard
#                     deviation that REML estimates for the 19 studies;
#   normal, SD 0.31   the same mean, and the posterior mean of the standard
#                     deviation of g in raudenbush_fit;
#   fitted g          the posterior mean of g in raudenbush_fit: a grid point
#                     drawn with its probability, moved uniformly within its
#                     cell of the grid.
#
# Each data set is analyzed in four ways:
#
#   bayesEfron                     bayes_efron_fit() with its defaults;
#   bayesEfron, adapt_delta = 0.99 the same with smaller steps of the sampler;
#   normal                         the normal model with its mean and standard
#                                  deviation held at their REML estimates for
#                                  that data set;
#   known g                        the posterior distribution under the
#                                  distribution that the effects were drawn
#                                  from, which no analysis of real data can
#                                  use.
#
# The defaults are 200 replications for each distribution and 12 worker
# processes. With a third argument the results of every data set are also
# written to that file. The seeds are fixed, so the table is the same when
# the script is run again with the same versions of CmdStan and the
# compiler.

args <- commandArgs(trailingOnly = TRUE)
analytic_only <- length(args) > 0L && identical(args[[1L]], "--analytic-only")
if (analytic_only) {
  rows_file <- if (length(args) >= 2L) args[[2L]] else
    file.path("data-raw", "benchmark", "results-v040", "few-studies-rows.csv")
  table_file <- if (length(args) >= 3L) args[[3L]] else
    file.path("inst", "verification", "few-studies.csv")
} else {
  replications <- if (length(args) >= 1L) as.integer(args[[1L]]) else 200L
  workers <- if (length(args) >= 2L) as.integer(args[[2L]]) else 12L
  rows_file <- if (length(args) >= 3L) args[[3L]] else NULL
  table_file <- file.path("inst", "verification", "few-studies.csv")
}
suppressMessages(devtools::load_all(".", quiet = TRUE))
source(file.path("data-raw", "benchmark", "histogram-posterior.R"))

sigma <- raudenbush1985$sei
K <- length(sigma)
reml <- .bef_reml_normal(raudenbush1985$yi, sigma)
# The posterior mean of g of the stored fit, on its grid.
fitted_g <- g_summary(raudenbush_fit, type = "mass")
grid <- fitted_g$value
mass <- fitted_g$mean / sum(fitted_g$mean)
spacing <- grid[2L] - grid[1L]
fitted_mean <- sum(mass * grid)
observed <- c(own = diagnose(raudenbush_fit)$grid$sd_g[["mean"]], reml = reml$sd)
truths <- data.frame(
  scenario = c("normal, SD 0.14", "normal, SD 0.31", "fitted g"),
  shape = c("normal", "normal", "fitted"),
  mean = c(reml$mean, reml$mean, fitted_mean),
  # A cell of the grid adds the variance of a uniform distribution on it.
  sd = c(reml$sd, observed[["own"]],
         sqrt(sum(mass * (grid - fitted_mean)^2) + spacing^2 / 12)),
  stringsAsFactors = FALSE
)
stopifnot(sprintf("normal, SD %.2f", truths$sd[1:2]) == truths$scenario[1:2])
level <- 0.9
z <- stats::qnorm(1 - (1 - level) / 2)
draw_effects <- function(truth) {
  if (identical(truth$shape, "normal")) {
    return(stats::rnorm(K, truth$mean, truth$sd))
  }
  sample(grid, K, replace = TRUE, prob = mass) + stats::runif(K, -spacing / 2, spacing / 2)
}

# The posterior means and 90% intervals of a normal model with a given mean
# and standard deviation of the effects.
normal_sites <- function(y, mu, tau) {
  weight <- tau^2 / (tau^2 + sigma^2)
  mean <- mu + weight * (y - mu)
  sd <- sqrt(weight) * sigma
  data.frame(mean = mean, lower = mean - z * sd, upper = mean + z * sd)
}
# The same under the fitted g, whose probability at a grid point is spread
# evenly over the cell of that point, as the effects are drawn.
fitted_sites <- function(y) {
  histogram_sites(y, sigma, grid, mass, level)
}
measures <- function(sites, theta, y) {
  largest <- which.max(y)
  inside <- theta >= sites$lower & theta <= sites$upper
  data.frame(
    inclusion = mean(inside), length = mean(sites$upper - sites$lower),
    rmse = sqrt(mean((sites$mean - theta)^2)),
    largest_inside = inside[largest],
    largest_error = sites$mean[largest] - theta[largest]
  )
}

# A fit of this package with the given arguments: the measures, and what the
# sampler and the grid say about it.
own_fit <- function(method, y, theta, truth, seed, ...) {
  warned <- character()
  fit <- withCallingHandlers(
    suppressMessages(bayes_efron_fit(y, sigma, seed = seed, parallel_chains = 1L, ...)),
    warning = function(w) {
      warned <<- c(warned, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }
  )
  d <- diagnose(fit)
  sd_g <- d$grid$sd_g
  cbind(
    method = method, measures(site_summary(fit, level = level), theta, y),
    sd = sd_g[["mean"]], sd_inside = truth$sd >= sd_g[["q5"]] & truth$sd <= sd_g[["q95"]],
    max_rhat = max(d$rhat, na.rm = TRUE), divergences = sum(d$divergences),
    end_warning = "bef_warn_grid_ends" %in% warned,
    width = diff(range(fit$metadata$data_list$grid)),
    stringsAsFactors = FALSE
  )
}
other_fit <- function(method, y, theta, sites, tau) {
  cbind(
    method = method, measures(sites, theta, y),
    sd = tau, sd_inside = NA, max_rhat = NA_real_, divergences = NA_real_,
    end_warning = NA, width = NA_real_, stringsAsFactors = FALSE
  )
}

one <- function(scenario, r) {
  truth <- truths[scenario, ]
  seed <- 20261200L + 1000L * scenario + r
  set.seed(seed)
  theta <- draw_effects(truth)
  y <- stats::rnorm(K, theta, sigma)
  estimate <- .bef_reml_normal(y, sigma)
  known <- if (identical(truth$shape, "normal")) {
    normal_sites(y, truth$mean, truth$sd)
  } else {
    fitted_sites(y)
  }
  rows <- rbind(
    own_fit("bayesEfron", y, theta, truth, seed),
    own_fit("bayesEfron, adapt_delta = 0.99", y, theta, truth, seed, adapt_delta = 0.99),
    other_fit("normal", y, theta, normal_sites(y, estimate$mean, estimate$sd), estimate$sd),
    other_fit("known g", y, theta, known, truth$sd)
  )
  cbind(scenario = truth$scenario, true_sd = truth$sd, replication = r, rows,
        stringsAsFactors = FALSE)
}

# Replace selected CSV records without rewriting any other line. The
# verification files have one record per line and no embedded newlines.
replace_csv_rows <- function(path, rows, selected) {
  lines <- readLines(path, warn = FALSE)
  stopifnot(length(lines) == nrow(rows) + 1L)
  for (i in which(selected)) {
    connection <- textConnection("record", "w", local = TRUE)
    utils::write.csv(rows[i, , drop = FALSE], connection, row.names = FALSE)
    close(connection)
    stopifnot(length(record) == 2L)
    lines[i + 1L] <- record[2L]
  }
  writeLines(lines, path, useBytes = TRUE)
}

if (analytic_only) {
  rows <- utils::read.csv(rows_file, stringsAsFactors = FALSE)
  rows_before <- rows
  known <- rows$scenario == "fitted g" & rows$method == "known g"
  stopifnot(any(known), !anyDuplicated(rows$replication[known]))
  truth <- truths[truths$shape == "fitted", ]
  stopifnot(all(abs(rows$true_sd[known] - truth$sd) < 1e-12))
  for (i in which(known)) {
    # The seed of this data set in one(): the third distribution.
    set.seed(20261200L + 3000L + rows$replication[i])
    theta <- draw_effects(truth)
    y <- stats::rnorm(K, theta, sigma)
    reference <- measures(fitted_sites(y), theta, y)
    rows[i, names(reference)] <- reference
  }
  stopifnot(identical(rows[!known, ], rows_before[!known, ]))
  message(sprintf("Computed %d known-g rows again; %d other rows are as they were.",
                  sum(known), sum(!known)))
} else {
  # One fit in this process first, so that the Stan program is compiled
  # before the workers start.
  invisible(one(1L, 0L))
  jobs <- expand.grid(r = seq_len(replications), scenario = seq_len(nrow(truths)))
  started <- Sys.time()
  results <- parallel::mclapply(
    seq_len(nrow(jobs)),
    function(i) tryCatch(one(jobs$scenario[i], jobs$r[i]), error = function(e) conditionMessage(e)),
    mc.cores = workers, mc.preschedule = FALSE, mc.silent = TRUE
  )
  failed <- vapply(results, is.character, logical(1))
  if (any(failed)) {
    stop(sprintf("%d fits failed; the first: %s", sum(failed), results[[which(failed)[1L]]]))
  }
  rows <- do.call(rbind, results)
  message(sprintf("%d data sets in %.1f minutes", nrow(jobs),
                  as.numeric(difftime(Sys.time(), started, units = "mins"))))
}

se <- function(x) stats::sd(x) / sqrt(length(x))
methods <- c("bayesEfron", "bayesEfron, adapt_delta = 0.99", "normal", "known g")
summary_of <- function(x) {
  diverged <- x$divergences > 0
  part <- function(keep) if (all(is.na(keep)) || !any(keep)) NA_real_ else mean(x$inclusion[keep])
  data.frame(
    scenario = x$scenario[1L], true_sd = x$true_sd[1L], method = x$method[1L],
    replications = nrow(x),
    # The estimate of the standard deviation of the effects: the posterior
    # mean of the standard deviation of g, or the REML estimate.
    sd_mean = mean(x$sd), sd_se = se(x$sd),
    sd_q05 = stats::quantile(x$sd, 0.05, names = FALSE),
    sd_q50 = stats::quantile(x$sd, 0.5, names = FALSE),
    sd_q95 = stats::quantile(x$sd, 0.95, names = FALSE),
    sd_zero = mean(x$sd == 0),
    sd_inclusion = if (all(is.na(x$sd_inside))) NA_real_ else mean(x$sd_inside),
    inclusion = mean(x$inclusion), inclusion_se = se(x$inclusion),
    length = mean(x$length), length_se = se(x$length),
    rmse = mean(x$rmse), rmse_se = se(x$rmse),
    # The study with the largest estimate.
    largest_inclusion = mean(x$largest_inside),
    largest_error = mean(x$largest_error), largest_error_se = se(x$largest_error),
    rhat_above_1.01 = sum(x$max_rhat > 1.01), with_divergences = sum(diverged),
    with_many_divergences = sum(x$divergences > 10),
    divergences_median = stats::median(x$divergences), divergences_max = max(x$divergences),
    inclusion_with_divergences = part(diverged), inclusion_without_divergences = part(!diverged),
    end_warnings = sum(x$end_warning), width_mean = mean(x$width),
    stringsAsFactors = FALSE
  )
}
groups <- split(rows, list(rows$method, rows$scenario), drop = TRUE)
table <- do.call(rbind, lapply(groups, summary_of))

# How often the two estimates of the standard deviation come out as they do
# for the 19 studies themselves: this package's as large, REML's as small.
wide <- merge(rows[rows$method == "bayesEfron", c("scenario", "replication", "sd")],
              rows[rows$method == "normal", c("scenario", "replication", "sd")],
              by = c("scenario", "replication"), suffixes = c("_own", "_reml"))
as_observed <- do.call(rbind, lapply(split(wide, wide$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1L],
    observed_own_sd = observed[["own"]], observed_reml_sd = observed[["reml"]],
    own_as_large = mean(x$sd_own >= observed[["own"]]),
    reml_as_small = mean(x$sd_reml <= observed[["reml"]]),
    both_as_observed = mean(x$sd_own >= observed[["own"]] & x$sd_reml <= observed[["reml"]]),
    stringsAsFactors = FALSE
  )
}))
table <- merge(table, as_observed, by = "scenario", sort = FALSE)

# The two fits of this package on the same data sets: what the smaller
# steps change, as the mean difference from the fit with the defaults.
paired <- merge(rows[rows$method == methods[1L], ], rows[rows$method == methods[2L], ],
                by = c("scenario", "replication"), suffixes = c("", "_99"))
changes <- do.call(rbind, lapply(split(paired, paired$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1L], method = methods[2L],
    sd_change = mean(x$sd_99 - x$sd), sd_change_se = se(x$sd_99 - x$sd),
    inclusion_change = mean(x$inclusion_99 - x$inclusion),
    inclusion_change_se = se(x$inclusion_99 - x$inclusion),
    length_change = mean(x$length_99 - x$length), length_change_se = se(x$length_99 - x$length),
    stringsAsFactors = FALSE
  )
}))
table <- merge(table, changes, by = c("scenario", "method"), all.x = TRUE, sort = FALSE)
table <- table[order(match(table$scenario, truths$scenario), match(table$method, methods)), ]
rownames(table) <- NULL
if (analytic_only) {
  table_before <- utils::read.csv(table_file, stringsAsFactors = FALSE)
  known_row <- table_before$scenario == "fitted g" & table_before$method == "known g"
  computed <- table[table$scenario == "fitted g" & table$method == "known g", ]
  stopifnot(sum(known_row) == 1L, nrow(computed) == 1L,
            table_before$replications[known_row] == sum(known))
  fields <- c("inclusion", "inclusion_se", "length", "length_se", "rmse", "rmse_se",
              "largest_inclusion", "largest_error", "largest_error_se")
  table <- table_before
  table[known_row, fields] <- computed[, fields]
  replace_csv_rows(table_file, table, known_row)
} else {
  table$generated_on <- format(Sys.Date())
  table$cmdstan_version <- raudenbush_fit$metadata$cmdstan_version
  utils::write.csv(table, table_file, row.names = FALSE)
}

options(width = 200)
print(table[, c("scenario", "method", "replications", "sd_mean", "sd_q05", "sd_q95", "sd_zero",
                "sd_inclusion", "inclusion", "inclusion_se", "length", "rmse")], digits = 3)
print(table[, c("method", "largest_inclusion", "largest_error", "rhat_above_1.01",
                "with_divergences", "with_many_divergences", "divergences_max",
                "inclusion_with_divergences", "inclusion_without_divergences", "end_warnings")],
      digits = 3)
print(table[, c("scenario", "method", "true_sd", "own_as_large", "reml_as_small",
                "both_as_observed", "sd_change", "inclusion_change", "length_change")],
      digits = 3)
if (!is.null(rows_file)) {
  if (analytic_only) {
    replace_csv_rows(rows_file, rows, known)
  } else {
    utils::write.csv(rows, rows_file, row.names = FALSE)
  }
}
