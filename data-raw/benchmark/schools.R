# Builds vignettes/articles/class-size-samples.csv: what a fit of this
# package and the normal random-effects model report for 79 schools when the
# truth is known. The article "Small classes in 79 schools" reads the table.
# Needs CmdStan and the stored fit and the evidence file of the article
# (data-raw/class-size.R). Run from the package root:
#
#   Rscript data-raw/benchmark/schools.R [replications] [workers] [rows-file]
#   Rscript data-raw/benchmark/schools.R --summary-only [rows-file]
#
# The second form makes no fits. It computes the table again from the
# results of every data set in the rows file, by default
# data-raw/benchmark/results-v040/schools-rows.csv, which the first form
# wrote.
#
# Every data set has the 79 standard errors of `star_schools`, and an
# estimate is its true effect plus normal noise with its standard error. The
# true effects come from one of four distributions:
#
#   no variation      every school has the mean that REML estimates for the
#                     79 schools;
#   normal, SD 0.14   a normal distribution with that mean and the standard
#                     deviation that REML estimates;
#   normal, SD 0.20   the same mean, and the posterior mean of the standard
#                     deviation of g in the stored fit of the article;
#   normal, SD 0.35   the same mean, and the upper limit of the 90% interval
#                     of that standard deviation.
#
# Each data set is analyzed in four ways:
#
#   bayesEfron             bayes_efron_fit() with the sampler settings of the
#                          stored fit: adapt_delta = 0.99, max_treedepth = 12,
#                          2,000 warmup and 6,000 sampling iterations. Every
#                          draw is kept, the chains run one after another,
#                          and each data set has its own seed;
#   bayesEfron, defaults   bayes_efron_fit() with its defaults;
#   normal                 the normal model with its mean and standard
#                          deviation held at their REML estimates for that
#                          data set;
#   known g                the posterior distribution under the distribution
#                          that the effects were drawn from, which no
#                          analysis of real data can use. Without variation
#                          it is the common effect itself, and those rows are
#                          left out.
#
# The defaults are 200 replications for each distribution and 12 worker
# processes. With a third argument the results of every data set are also
# written to that file. The seeds are fixed, so the table is the same when
# the script is run again with the same versions of CmdStan and the
# compiler.

args <- commandArgs(trailingOnly = TRUE)
summary_only <- length(args) > 0L && identical(args[[1L]], "--summary-only")
if (summary_only) {
  rows_file <- if (length(args) >= 2L) args[[2L]] else
    file.path("data-raw", "benchmark", "results-v040", "schools-rows.csv")
} else {
  replications <- if (length(args) >= 1L) as.integer(args[[1L]]) else 200L
  workers <- if (length(args) >= 2L) as.integer(args[[2L]]) else 12L
  rows_file <- if (length(args) >= 3L) args[[3L]] else NULL
}
table_file <- file.path("vignettes", "articles", "class-size-samples.csv")
suppressMessages(devtools::load_all(".", quiet = TRUE))

stored <- readRDS(file.path("vignettes", "articles", "class-size-fit.rds"))
stopifnot(identical(stored$metadata$data_list$theta_hat, star_schools$yi),
          identical(stored$metadata$data_list$sigma, star_schools$sei))
sigma <- star_schools$sei
K <- length(sigma)
reml <- .bef_reml_normal(star_schools$yi, sigma)
# The posterior mean of the standard deviation of g: in the stored fit, which
# keeps 800 draws, and in all the 24,000 draws of its run.
evidence <- readRDS(file.path("vignettes", "articles", "class-size-evidence.rds"))
all_draws <- evidence$runs$sd_g[evidence$runs$setting == "stored" & evidence$runs$seed == 1985L]
stopifnot(length(all_draws) == 1L)
observed <- c(own = diagnose(stored)$grid$sd_g[["mean"]], reml = reml$sd, own_all_draws = all_draws)
truths <- data.frame(
  scenario = c("no variation", "normal, SD 0.14", "normal, SD 0.20", "normal, SD 0.35"),
  mean = reml$mean,
  sd = c(0, reml$sd, observed[["own"]], diagnose(stored)$grid$sd_g[["q95"]]),
  stringsAsFactors = FALSE
)
stopifnot(sprintf("normal, SD %.2f", truths$sd[2:4]) == truths$scenario[2:4])
level <- 0.9
z <- stats::qnorm(1 - (1 - level) / 2)
stored_settings <- list(adapt_delta = 0.99, max_treedepth = 12L,
                        iter_warmup = 2000L, iter_sampling = 6000L)

# The posterior means and 90% intervals of a normal model with a given mean
# and standard deviation of the effects.
normal_sites <- function(y, mu, tau) {
  weight <- tau^2 / (tau^2 + sigma^2)
  mean <- mu + weight * (y - mu)
  sd <- sqrt(weight) * sigma
  data.frame(mean = mean, lower = mean - z * sd, upper = mean + z * sd)
}
measures <- function(sites, theta, y) {
  largest <- which.max(y)
  inside <- theta >= sites$lower & theta <= sites$upper
  data.frame(
    inclusion = mean(inside), length = mean(sites$upper - sites$lower),
    rmse = sqrt(mean((sites$mean - theta)^2)),
    largest_inside = inside[largest],
    largest_error = sites$mean[largest] - theta[largest],
    excluding_zero = sum(sites$lower > 0 | sites$upper < 0)
  )
}

# A fit of this package with the given arguments: the measures, and what the
# sampler and the grid say about it.
own_fit <- function(method, y, theta, truth, seed, settings) {
  warned <- character()
  fit <- withCallingHandlers(
    suppressMessages(do.call(
      bayes_efron_fit, c(list(y, sigma, seed = seed, parallel_chains = 1L), settings)
    )),
    warning = function(w) {
      warned <<- c(warned, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }
  )
  d <- diagnose(fit)
  sd_g <- d$grid$sd_g
  g <- confint(fit, type = "g", level = level)
  rownames(g) <- g$parameter
  cbind(
    method = method, measures(site_summary(fit, level = level), theta, y),
    mean = g["mean_g", "mean"],
    mean_inside = truth$mean >= g["mean_g", "lower"] & truth$mean <= g["mean_g", "upper"],
    sd = sd_g[["mean"]], sd_lower = sd_g[["q5"]], sd_upper = sd_g[["q95"]],
    sd_inside = truth$sd >= sd_g[["q5"]] & truth$sd <= sd_g[["q95"]],
    max_rhat = max(d$rhat, na.rm = TRUE), min_ess_bulk = min(d$ess_bulk, na.rm = TRUE),
    ebfmi_min = min(d$ebfmi), divergences = sum(d$divergences),
    treedepth_hits = sum(d$max_treedepth),
    end_warning = "bef_warn_grid_ends" %in% warned,
    stringsAsFactors = FALSE
  )
}
other_fit <- function(method, y, theta, sites, mu, tau) {
  cbind(
    method = method, measures(sites, theta, y),
    mean = mu, mean_inside = NA, sd = tau, sd_lower = NA_real_, sd_upper = NA_real_,
    sd_inside = NA, max_rhat = NA_real_, min_ess_bulk = NA_real_, ebfmi_min = NA_real_,
    divergences = NA_real_, treedepth_hits = NA_real_, end_warning = NA,
    stringsAsFactors = FALSE
  )
}

one <- function(scenario, r) {
  truth <- truths[scenario, ]
  seed <- 20261300L + 1000L * scenario + r
  set.seed(seed)
  theta <- stats::rnorm(K, truth$mean, truth$sd)
  y <- stats::rnorm(K, theta, sigma)
  estimate <- .bef_reml_normal(y, sigma)
  rows <- rbind(
    own_fit("bayesEfron", y, theta, truth, seed, stored_settings),
    own_fit("bayesEfron, defaults", y, theta, truth, seed, list()),
    other_fit("normal", y, theta, normal_sites(y, estimate$mean, estimate$sd),
              estimate$mean, estimate$sd)
  )
  if (truth$sd > 0) {
    rows <- rbind(rows, other_fit("known g", y, theta, normal_sites(y, truth$mean, truth$sd),
                                  truth$mean, truth$sd))
  }
  cbind(scenario = truth$scenario, true_mean = truth$mean, true_sd = truth$sd,
        replication = r, rows, stringsAsFactors = FALSE)
}

if (summary_only) {
  rows <- utils::read.csv(rows_file, stringsAsFactors = FALSE)
  stopifnot(identical(unique(rows$scenario), truths$scenario),
            all(abs(rows$true_sd - truths$sd[match(rows$scenario, truths$scenario)]) < 1e-12))
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
  rows$generated_on <- format(Sys.Date())
  rows$cmdstan_version <- as.character(cmdstanr::cmdstan_version())
  message(sprintf("%d data sets in %.1f minutes", nrow(jobs),
                  as.numeric(difftime(Sys.time(), started, units = "mins"))))
}

se <- function(x) stats::sd(x) / sqrt(length(x))
share <- function(x) if (all(is.na(x))) NA_real_ else mean(x)
count <- function(x) if (all(is.na(x))) NA_integer_ else sum(x)
methods <- c("bayesEfron", "bayesEfron, defaults", "normal", "known g")
summary_of <- function(x) {
  low <- x$ebfmi_min < 0.3
  part <- function(keep) if (all(is.na(keep)) || !any(keep)) NA_real_ else mean(x$inclusion[keep])
  # The data sets whose fit reports a standard deviation no larger than that
  # of the schools.
  as_small <- if (x$method[1L] %in% methods[1:2]) x$sd <= observed[["own"]] else NA
  data.frame(
    scenario = x$scenario[1L], true_mean = x$true_mean[1L], true_sd = x$true_sd[1L],
    method = x$method[1L], replications = nrow(x),
    # The estimate of the mean of the effects: the posterior mean of the mean
    # of g, or the REML estimate.
    mean_mean = mean(x$mean), mean_inclusion = share(x$mean_inside),
    # The estimate of the standard deviation of the effects: the posterior
    # mean of the standard deviation of g, or the REML estimate.
    sd_mean = mean(x$sd), sd_se = se(x$sd),
    sd_q05 = stats::quantile(x$sd, 0.05, names = FALSE),
    sd_q50 = stats::quantile(x$sd, 0.5, names = FALSE),
    sd_q95 = stats::quantile(x$sd, 0.95, names = FALSE),
    sd_zero = mean(x$sd == 0),
    sd_lower_mean = mean(x$sd_lower), sd_upper_mean = mean(x$sd_upper),
    sd_inclusion = share(x$sd_inside),
    inclusion = mean(x$inclusion), inclusion_se = se(x$inclusion),
    sd_as_small = count(as_small), inclusion_sd_as_small = part(as_small),
    length = mean(x$length), length_se = se(x$length),
    rmse = mean(x$rmse), rmse_se = se(x$rmse),
    # The school with the largest estimate.
    largest_inclusion = mean(x$largest_inside),
    largest_error = mean(x$largest_error), largest_error_se = se(x$largest_error),
    excluding_zero_mean = mean(x$excluding_zero),
    excluding_zero_q05 = stats::quantile(x$excluding_zero, 0.05, names = FALSE),
    excluding_zero_q95 = stats::quantile(x$excluding_zero, 0.95, names = FALSE),
    ebfmi_below_0.3 = count(low), ebfmi_median = stats::median(x$ebfmi_min),
    inclusion_ebfmi_below = part(low), inclusion_ebfmi_above = part(!low),
    # The E-BFMI of a fit against its estimate of the standard deviation.
    ebfmi_sd_correlation = if (all(is.na(low))) NA_real_ else stats::cor(x$ebfmi_min, x$sd),
    rhat_above_1.01 = count(x$max_rhat > 1.01), ess_below_400 = count(x$min_ess_bulk < 400),
    with_divergences = count(x$divergences > 0),
    with_treedepth_hits = count(x$treedepth_hits > 0),
    end_warnings = count(x$end_warning),
    stringsAsFactors = FALSE
  )
}
groups <- split(rows, list(rows$method, rows$scenario), drop = TRUE)
table <- do.call(rbind, lapply(groups, summary_of))

# How often the two estimates of the standard deviation come out at least as
# large as they do for the 79 schools themselves.
wide <- merge(rows[rows$method == "bayesEfron", c("scenario", "replication", "sd")],
              rows[rows$method == "normal", c("scenario", "replication", "sd")],
              by = c("scenario", "replication"), suffixes = c("_own", "_reml"))
as_observed <- do.call(rbind, lapply(split(wide, wide$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1L],
    observed_own_sd = observed[["own"]], observed_reml_sd = observed[["reml"]],
    own_as_large = mean(x$sd_own >= observed[["own"]]),
    own_as_large_all_draws = mean(x$sd_own >= observed[["own_all_draws"]]),
    reml_as_large = mean(x$sd_reml >= observed[["reml"]]),
    stringsAsFactors = FALSE
  )
}))
table <- merge(table, as_observed, by = "scenario", sort = FALSE)

# The two fits of this package on the same data sets: what the defaults
# change, as the mean difference from the fit with the settings of the
# stored fit.
paired <- merge(rows[rows$method == methods[1L], ], rows[rows$method == methods[2L], ],
                by = c("scenario", "replication"), suffixes = c("", "_default"))
changes <- do.call(rbind, lapply(split(paired, paired$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1L], method = methods[2L],
    sd_change = mean(x$sd_default - x$sd), sd_change_largest = max(abs(x$sd_default - x$sd)),
    inclusion_change = mean(x$inclusion_default - x$inclusion),
    inclusion_change_largest = max(abs(x$inclusion_default - x$inclusion)),
    length_change = mean(x$length_default - x$length),
    stringsAsFactors = FALSE
  )
}))
table <- merge(table, changes, by = c("scenario", "method"), all.x = TRUE, sort = FALSE)
table <- table[order(match(table$scenario, truths$scenario), match(table$method, methods)), ]
rownames(table) <- NULL
# The date and the version of CmdStan of the fits, which the rows carry.
stopifnot(length(unique(rows$generated_on)) == 1L, length(unique(rows$cmdstan_version)) == 1L)
table$generated_on <- rows$generated_on[1L]
table$cmdstan_version <- rows$cmdstan_version[1L]
utils::write.csv(table, table_file, row.names = FALSE)

options(width = 220)
print(table[, c("scenario", "method", "replications", "mean_mean", "mean_inclusion", "sd_mean",
                "sd_q05", "sd_q95", "sd_zero", "sd_inclusion", "inclusion", "inclusion_se",
                "length", "rmse")], digits = 3)
print(table[, c("scenario", "method", "largest_inclusion", "largest_error", "excluding_zero_mean",
                "excluding_zero_q05", "excluding_zero_q95", "ebfmi_below_0.3", "ebfmi_median",
                "rhat_above_1.01", "ess_below_400", "with_divergences", "with_treedepth_hits",
                "end_warnings")], digits = 3)
print(table[, c("scenario", "method", "observed_own_sd", "observed_reml_sd", "own_as_large",
                "reml_as_large", "sd_lower_mean", "sd_upper_mean", "sd_change",
                "sd_change_largest", "inclusion_change", "length_change")], digits = 3)
if (!summary_only && !is.null(rows_file)) utils::write.csv(rows, rows_file, row.names = FALSE)
