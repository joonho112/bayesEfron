# Builds inst/examples/ebnm_comparison.rds: the stored fits of this package
# beside empirical Bayes fits of the same data by the ebnm package, and
# beside Efron's estimator where the standard errors are equal. The
# documentation reads the file, so that it is built without ebnm. Needs
# ebnm, deconvolveR and CmdStan. Run from the package root after
# data-raw/example-fits.R and the panel runs named below:
#
#   Rscript data-raw/ebnm-comparison.R [audit-dir]
#
# The methods, their settings and the measures are those of
# data-raw/ebnm-protocol.md, which was fixed before this script
# was run. The parts:
#
#   krw       the 97 callback gaps: agreement between the methods, the ten
#             firms with the largest posterior means under each, and the
#             distribution functions of the estimates of g;
#   twin      the 1,500 simulated sites, whose true effects are known:
#             error, inclusion, length, and distance to the truth;
#   equal_se  the same true effects with standard errors of 1: Efron's
#             estimator on the grid and with the number of spline functions
#             of the fit of this package; and the same over 20 such data
#             sets (data-raw/benchmark/equal-se-run.R);
#   panels    20 independent panels of 1,500 sites and 100 of 100 sites:
#             error, inclusion and length with their Monte Carlo errors,
#             and the convergence of the fits of the large panels. These
#             two tables are also written to inst/verification/.
#
# The fits of the panels are read from data-raw/benchmark/results-v040/,
# which data-raw/benchmark/panels-run.R wrote (arms "default", "M8" and
# "v030" for 1,500 sites, "default" for 100). Two fits are made here: the
# simulated sites on the range around the true effects, and the data with
# equal standard errors. With an audit directory the tables are also
# written there as CSV files. The times are those of this machine; all
# other numbers are the same when the script is run again.

args <- commandArgs(trailingOnly = TRUE)
audit_dir <- if (length(args) >= 1L) args[[1L]] else NULL
suppressMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("ebnm", quietly = TRUE), requireNamespace("deconvolveR", quietly = TRUE))

probs <- c(0.05, 0.95)
n_sim <- 4000L
everything <- ebnm::ebnm_output_all()

# ---- the empirical Bayes methods ---------------------------------------------
methods <- list(
  "No shrinkage" = function(x, s) ebnm::ebnm_flat(x, s, output = everything),
  "Normal" = function(x, s) ebnm::ebnm_normal(x, s, mode = "estimate", output = everything),
  "Scale mixture of normals" = function(x, s) {
    ebnm::ebnm_normal_scale_mixture(x, s, mode = "estimate", output = everything)
  },
  "ash" = function(x, s) {
    ebnm::ebnm_ash(x, s, mode = "estimate", mixcompdist = "normal", prior = "uniform",
                   output = everything)
  },
  "Unimodal" = function(x, s) ebnm::ebnm_unimodal(x, s, mode = "estimate", output = everything),
  "NPMLE" = function(x, s) ebnm::ebnm_npmle(x, s, output = everything)
)
calls <- c(
  "No shrinkage" = "ebnm_flat(x, s)",
  "Normal" = "ebnm_normal(x, s, mode = \"estimate\")",
  "Scale mixture of normals" = "ebnm_normal_scale_mixture(x, s, mode = \"estimate\")",
  "ash" = "ebnm_ash(x, s, mode = \"estimate\", mixcompdist = \"normal\", prior = \"uniform\")",
  "Unimodal" = "ebnm_unimodal(x, s, mode = \"estimate\")",
  "NPMLE" = "ebnm_npmle(x, s)",
  "Known g" = "ebnm_npmle(x, s, g_init = <true g>, fix_g = TRUE)",
  "g of bayesEfron, plugged in" = "ebnm_npmle(x, s, g_init = as_normalmix(fit), fix_g = TRUE)",
  "bayesEfron" = "bayes_efron_fit()",
  "bayesEfron, range from the true effects" =
    "bayes_efron_fit(grid_method = \"paper_simulation\", theta_true = , M = 8)",
  "Efron's estimator" = "deconvolveR::deconv(tau, P, Q, y = 1, family = \"Binomial\", c0 = 1)"
)

# ebnm says for every fit with an estimated mode that it returns no local
# false sign rates; nothing here uses them.
quietly <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("local false sign rates", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# What is taken from an ebnm fit: posterior means and standard deviations,
# the 5% and 95% sample quantiles of 4,000 draws from each posterior, the
# maximized log likelihood, the prior and the time.
fit_ebnm <- function(fun, x, s) {
  started <- Sys.time()
  fit <- quietly(fun(x, s))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  limits <- stats::quantile(fit, probs = probs, nsim = n_sim, seed = 1L)
  list(mean = unname(stats::coef(fit)), sd = stats::fitted(fit)$sd,
       lower = unname(limits[, 1L]), upper = unname(limits[, 2L]),
       loglik = as.numeric(stats::logLik(fit)), g = fit$fitted_g, seconds = seconds)
}

# The distribution function of a prior of ebnm at the values `x`: a mixture
# of normals or of point masses, or of uniforms. The flat prior has none.
cdf_of <- function(g, x) {
  if (inherits(g, "normalmix")) {
    if (any(!is.finite(g$sd))) {
      return(rep(NA_real_, length(x)))
    }
    smooth <- g$sd > 0
    return(vapply(x, function(v) {
      sum(g$pi[smooth] * stats::pnorm(v, g$mean[smooth], g$sd[smooth])) +
        sum(g$pi[!smooth] * (v >= g$mean[!smooth]))
    }, numeric(1L)))
  }
  if (inherits(g, "unimix")) {
    wide <- g$b > g$a
    return(vapply(x, function(v) {
      sum(g$pi[wide] * stats::punif(v, g$a[wide], g$b[wide])) +
        sum(g$pi[!wide] * (v >= g$a[!wide]))
    }, numeric(1L)))
  }
  rep(NA_real_, length(x))
}

# The posterior means, standard deviations and 90% intervals of the sites
# under one g on a grid, computed as the package computes them for a fit.
under_one_g <- function(grid, g, x, s, tau_c = 0) {
  engine <- .bef_engine_from_parts(log_g = matrix(log(g / sum(g)), nrow = 1L), grid = grid,
                                   theta_hat = x, sigma = s, tau_c = tau_c)
  weights <- .bef_site_weights(engine)
  limits <- .bef_site_quantile(engine, probs, weights = weights)
  moments <- .bef_site_moments(engine)
  list(mean = moments$mean, sd = sqrt(moments$var), lower = limits[, 1L], upper = limits[, 2L],
       engine = engine)
}

# The estimate of g of a fit of this package as a distribution function,
# its point masses spread over their cells as the package reads them.
cdf_of_fit_g <- function(grid, g, x, tau_c = 0) {
  engine <- .bef_engine_from_parts(log_g = matrix(log(g / sum(g)), nrow = 1L), grid = grid,
                                   theta_hat = 0, sigma = 1, tau_c = tau_c)
  as.numeric((g / sum(g)) %*% .bef_point_mass_below(engine, x, "interpolated"))
}

# How well posterior means and intervals recover known effects.
recovery <- function(mean, lower, upper, truth, s) {
  third <- cut(rank(s, ties.method = "first"), 3L, labels = FALSE)
  rmse <- function(rows) sqrt(mean((mean[rows] - truth[rows])^2))
  inside <- truth >= lower & truth <= upper
  data.frame(
    rmse = rmse(TRUE), rmse_precise = rmse(third == 1L), rmse_middle = rmse(third == 2L),
    rmse_imprecise = rmse(third == 3L), inclusion = mean(inside),
    inclusion_se = sqrt(mean(inside) * (1 - mean(inside)) / length(truth)),
    mean_length = mean(upper - lower)
  )
}

# The Wasserstein distance between a distribution function on the values
# `x`, equally spaced, and the empirical distribution of `truth`.
w1 <- function(cdf, x, truth) {
  if (anyNA(cdf)) NA_real_ else sum(abs(cdf - stats::ecdf(truth)(x))) * (x[2L] - x[1L])
}

sampling_seconds <- function(fit) {
  pilot <- attr(fit$metadata, "grid", exact = TRUE)$pilot$seconds
  fit$metadata$runtime_seconds + if (is.null(pilot)) 0 else pilot
}

quiet_fit <- function(...) {
  suppressWarnings(suppressMessages(bayes_efron_fit(...)))
}

# ---- krw: the callback gaps --------------------------------------------------
message("krw")
x <- krw_discrimination$yi
s <- krw_discrimination$sei
firms <- krw_discrimination$firm_name
fit <- krw_fit
sites <- site_summary(fit)
stopifnot(identical(sites$label, firms))
fit_grid <- fit$metadata$data_list$grid
fit_g <- as_normalmix(fit)
rows <- lapply(methods, fit_ebnm, x = x, s = s)
plugged <- fit_ebnm(function(x, s) {
  ebnm::ebnm_npmle(x, s, g_init = fit_g, fix_g = TRUE, output = everything)
}, x, s)
one_g <- under_one_g(fit_grid, fit_g$pi, x, s)
stopifnot(isTRUE(all.equal(plugged$mean, one_g$mean, tolerance = 1e-8)))
plugged[c("lower", "upper")] <- one_g[c("lower", "upper")]
rows[["g of bayesEfron, plugged in"]] <- plugged
rows[["bayesEfron"]] <- list(mean = sites$mean, sd = sites$sd, lower = sites$lower,
                             upper = sites$upper, loglik = NA_real_,
                             seconds = sampling_seconds(fit))

krw_table <- do.call(rbind, lapply(names(rows), function(method) {
  row <- rows[[method]]
  data.frame(
    method = method,
    correlation = stats::cor(row$mean, sites$mean),
    largest_difference_in_se = max(abs(row$mean - sites$mean) / s),
    mean_sd = mean(row$sd), sd_ratio = mean(row$sd) / mean(sites$sd),
    mean_length = mean(row$upper - row$lower),
    length_ratio = mean(row$upper - row$lower) / mean(sites$upper - sites$lower),
    loglik = row$loglik, seconds = row$seconds, stringsAsFactors = FALSE
  )
}))
krw_firms <- do.call(rbind, lapply(names(rows), function(method) {
  top <- order(-rows[[method]]$mean)[1:10]
  data.frame(method = method, rank = 1:10, site = top, firm = firms[top],
             mean = rows[[method]]$mean[top], stringsAsFactors = FALSE)
}))
at <- seq(fit_grid[1L], fit_grid[length(fit_grid)], length.out = 401L)
band <- g_summary(fit, type = "cdf", at = at)
krw_cdf <- rbind(
  do.call(rbind, lapply(names(methods), function(method) {
    data.frame(method = method, x = at, cdf = cdf_of(rows[[method]]$g, at), lower = NA_real_,
               upper = NA_real_, stringsAsFactors = FALSE)
  })),
  data.frame(method = "bayesEfron", x = at, cdf = band$mean, lower = band$lower,
             upper = band$upper, stringsAsFactors = FALSE)
)
krw_cdf <- krw_cdf[!is.na(krw_cdf$cdf), ]

# ---- twin: the simulated sites -----------------------------------------------
message("twin")
x <- twin_towers$yi
s <- twin_towers$sei
truth <- twin_towers$theta_true
fit <- twin_towers_fit
sites <- site_summary(fit)
fit_grid <- fit$metadata$data_list$grid
fit_g <- as_normalmix(fit)
true_g <- structure(data.frame(pi = rep(1 / length(truth), length(truth)), mean = truth, sd = 0),
                    class = "normalmix")
rows <- lapply(methods, fit_ebnm, x = x, s = s)
rows[["Known g"]] <- fit_ebnm(function(x, s) {
  ebnm::ebnm_npmle(x, s, g_init = true_g, fix_g = TRUE, output = everything)
}, x, s)
plugged <- fit_ebnm(function(x, s) {
  ebnm::ebnm_npmle(x, s, g_init = fit_g, fix_g = TRUE, output = everything)
}, x, s)
one_g <- under_one_g(fit_grid, fit_g$pi, x, s)
stopifnot(isTRUE(all.equal(plugged$mean, one_g$mean, tolerance = 1e-8)))
plugged[c("lower", "upper")] <- one_g[c("lower", "upper")]
rows[["g of bayesEfron, plugged in"]] <- plugged
rows[["bayesEfron"]] <- list(mean = sites$mean, sd = sites$sd, lower = sites$lower,
                             upper = sites$upper, loglik = NA_real_,
                             seconds = sampling_seconds(fit))
oracle <- quiet_fit(twin_towers, grid_method = "paper_simulation", theta_true = truth, M = 8L,
                    seed = 20250627L, thin = 15L)
oracle_sites <- site_summary(oracle)
oracle_g <- as_normalmix(oracle)
rows[["bayesEfron, range from the true effects"]] <- list(
  mean = oracle_sites$mean, sd = oracle_sites$sd, lower = oracle_sites$lower,
  upper = oracle_sites$upper, loglik = NA_real_, seconds = sampling_seconds(oracle)
)

# The estimates of g as distribution functions on one fine support.
support <- seq(min(fit_grid[1L], oracle_g$mean[1L], x) - 1,
               max(fit_grid[length(fit_grid)], x) + 1, length.out = 40001L)
distribution <- function(method, on) {
  switch(
    method,
    "bayesEfron" = , "g of bayesEfron, plugged in" = cdf_of_fit_g(fit_grid, fit_g$pi, on),
    "bayesEfron, range from the true effects" = cdf_of_fit_g(oracle_g$mean, oracle_g$pi, on),
    cdf_of(rows[[method]]$g, on)
  )
}
# The log likelihood of the two fits of this package is that of their
# estimate of g taken as known.
rows[["bayesEfron, range from the true effects"]]$plug_in_loglik <- as.numeric(stats::logLik(
  ebnm::ebnm_npmle(x, s, g_init = oracle_g, fix_g = TRUE)
))
twin_table <- do.call(rbind, lapply(names(rows), function(method) {
  row <- rows[[method]]
  cbind(
    data.frame(method = method, stringsAsFactors = FALSE),
    recovery(row$mean, row$lower, row$upper, truth, s),
    data.frame(
      w1 = w1(distribution(method, support), support, truth),
      loglik = if (is.null(row$plug_in_loglik)) row$loglik else row$plug_in_loglik,
      seconds = row$seconds
    )
  )
}))
twin_table$loglik[twin_table$method == "bayesEfron"] <-
  twin_table$loglik[twin_table$method == "g of bayesEfron, plugged in"]
at <- seq(fit_grid[1L], fit_grid[length(fit_grid)], length.out = 401L)
band <- g_summary(fit, type = "cdf", at = at)
twin_cdf <- rbind(
  do.call(rbind, lapply(c(names(methods), "Known g"), function(method) {
    data.frame(method = method, x = at, cdf = cdf_of(rows[[method]]$g, at), lower = NA_real_,
               upper = NA_real_, stringsAsFactors = FALSE)
  })),
  data.frame(method = "bayesEfron", x = at, cdf = band$mean, lower = band$lower,
             upper = band$upper, stringsAsFactors = FALSE)
)
twin_cdf <- twin_cdf[!is.na(twin_cdf$cdf), ]

# ---- equal_se: Efron's estimator under its own conditions --------------------
message("equal_se")
K <- length(truth)
set.seed(20261101L)
z <- truth + stats::rnorm(K)
ones <- rep(1, K)
equal_fit <- quiet_fit(z, ones, M = 8L, seed = 20261102L, thin = 15L)
equal_sites <- site_summary(equal_fit)
tau <- equal_fit$metadata$data_list$grid
M <- equal_fit$metadata$data_list$M
# The basis that deconv() builds for the normal family, and the likelihood
# of every estimate at every grid point. deconv() called with the estimates
# alone would count them in 39 intervals first. With `P`, `Q` and `y = 1`,
# `family = "Binomial"` is the path of deconv() that takes one row of `P`
# for each observation; it does not change the likelihood, which is `P`.
Q <- scale(splines::ns(tau, M), center = TRUE, scale = FALSE)
Q <- apply(Q, 2L, function(w) w / sqrt(sum(w * w)))
P <- outer(z, tau, function(value, point) stats::dnorm(value - point))
started <- Sys.time()
efron <- deconvolveR::deconv(tau = tau, P = P, Q = Q, y = 1, family = "Binomial", c0 = 1)
efron_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
efron_g <- efron$stats[, "g"]
efron_rows <- under_one_g(tau, efron_g, z, ones)
equal_g <- as_normalmix(equal_fit)
equal_plug <- under_one_g(tau, equal_g$pi, z, ones)
support <- seq(min(tau[1L], z) - 1, max(tau[length(tau)], z) + 1, length.out = 40001L)
loglik_of <- function(g) sum(log(as.numeric(P %*% (g / sum(g)))))
equal_rows <- list(
  "Efron's estimator" = c(efron_rows, list(g = efron_g, seconds = efron_seconds)),
  "g of bayesEfron, plugged in" = c(equal_plug, list(g = equal_g$pi, seconds = NA_real_)),
  "bayesEfron" = list(mean = equal_sites$mean, sd = equal_sites$sd, lower = equal_sites$lower,
                      upper = equal_sites$upper, g = equal_g$pi,
                      seconds = sampling_seconds(equal_fit))
)
equal_table <- do.call(rbind, lapply(names(equal_rows), function(method) {
  row <- equal_rows[[method]]
  cbind(
    data.frame(method = method, stringsAsFactors = FALSE),
    recovery(row$mean, row$lower, row$upper, truth, ones),
    data.frame(mean_sd = mean(row$sd), w1 = w1(cdf_of_fit_g(tau, row$g, support), support, truth),
               loglik = loglik_of(row$g), seconds = row$seconds)
  )
}))
band <- g_summary(equal_fit, type = "density")
equal_g_table <- data.frame(
  x = tau, bayesEfron = band$mean, lower = band$lower, upper = band$upper,
  efron = efron_g / diff(tau[1:2]), efron_se = efron$stats[, "SE.g"] / diff(tau[1:2])
)
equal_agreement <- data.frame(
  grid_lower = tau[1L], grid_upper = tau[length(tau)], L = length(tau), M = M,
  largest_cdf_difference = max(abs(cumsum(efron_g / sum(efron_g)) - cumsum(equal_g$pi))),
  correlation_of_means = stats::cor(efron_rows$mean, equal_sites$mean),
  largest_difference_of_means = max(abs(efron_rows$mean - equal_sites$mean)),
  length_ratio = mean(equal_sites$upper - equal_sites$lower) /
    mean(efron_rows$upper - efron_rows$lower)
)

# ---- panels: repeated samples ------------------------------------------------
message("panels")
results <- file.path("data-raw", "benchmark", "results-v040")
over_panels <- function(values) c(mean = mean(values), se = stats::sd(values) / sqrt(length(values)))

# Equal standard errors on the true effects of the 20 large panels, by
# data-raw/benchmark/equal-se-run.R, for six and eight spline functions.
equal_runs <- do.call(rbind, lapply(c(6L, 8L), function(M) {
  utils::read.csv(file.path(results, sprintf("equal-se-M%d.csv", M)), stringsAsFactors = FALSE)
}))
equal_runs <- equal_runs[order(equal_runs$M, equal_runs$panel), ]
equal_order <- c("bayesEfron", "g of bayesEfron, plugged in", "Efron's estimator")
equal_repeated <- do.call(rbind, lapply(split(equal_runs, list(equal_runs$method, equal_runs$M)),
                                        function(rows) {
  fits <- rows[!duplicated(rows$panel), ]
  data.frame(
    M = rows$M[[1L]], method = rows$method[[1L]], panels = nrow(rows),
    rmse = over_panels(rows$rmse)[["mean"]], rmse_se = over_panels(rows$rmse)[["se"]],
    inclusion = over_panels(rows$inclusion)[["mean"]],
    inclusion_se = over_panels(rows$inclusion)[["se"]],
    mean_length = over_panels(rows$mean_length)[["mean"]],
    mean_length_se = over_panels(rows$mean_length)[["se"]],
    w1_population = over_panels(rows$w1_population)[["mean"]],
    w1_population_se = over_panels(rows$w1_population)[["se"]],
    rhat_max = max(fits$max_rhat_parameters), stringsAsFactors = FALSE
  )
}))
equal_repeated <- equal_repeated[order(equal_repeated$M,
                                       match(equal_repeated$method, equal_order)), ]
rownames(equal_repeated) <- NULL

panel_rows <- list()
for (K in c(100L, 1500L)) {
  panels <- readRDS(file.path("data-raw", "benchmark", sprintf("panels_K%d.rds", K)))$panels
  kept <- readRDS(file.path(results, sprintf("panels-K%d-default-fits.rds", K)))
  stopifnot(length(kept) == length(panels))
  for (p in seq_along(panels)) {
    panel <- panels[[p]]
    stored <- kept[[as.character(p)]]
    stopifnot(identical(stored$panel, panel$panel), nrow(stored$sites) == panel$K)
    x <- panel$theta_hat
    s <- panel$sigma
    one_g <- under_one_g(stored$grid, unname(stored$g), x, s)
    by_method <- list(
      "bayesEfron" = stored$sites,
      "g of bayesEfron, plugged in" = one_g,
      "Normal" = fit_ebnm(methods[["Normal"]], x, s),
      "NPMLE" = fit_ebnm(methods[["NPMLE"]], x, s)
    )
    for (method in names(by_method)) {
      row <- by_method[[method]]
      panel_rows[[length(panel_rows) + 1L]] <- cbind(
        data.frame(K = K, panel = p, method = method, stringsAsFactors = FALSE),
        recovery(row$mean, row$lower, row$upper, panel$theta_true, s)[
          c("rmse", "inclusion", "mean_length")]
      )
    }
  }
  message("  K = ", K, ": ", length(panels), " panels")
}
panel_rows <- do.call(rbind, panel_rows)
panel_table <- do.call(rbind, lapply(split(panel_rows, list(panel_rows$method, panel_rows$K)),
                                     function(rows) {
  data.frame(
    K = rows$K[[1L]], method = rows$method[[1L]], panels = nrow(rows),
    rmse = over_panels(rows$rmse)[["mean"]], rmse_se = over_panels(rows$rmse)[["se"]],
    inclusion = over_panels(rows$inclusion)[["mean"]],
    inclusion_se = over_panels(rows$inclusion)[["se"]],
    mean_length = over_panels(rows$mean_length)[["mean"]],
    mean_length_se = over_panels(rows$mean_length)[["se"]], stringsAsFactors = FALSE
  )
}))
order_of <- c("bayesEfron", "g of bayesEfron, plugged in", "Normal", "NPMLE")
panel_table <- panel_table[order(panel_table$K, match(panel_table$method, order_of)), ]
rownames(panel_table) <- NULL

# The convergence of the fits of the large panels, by setting.
arms <- c(default = "defaults", M8 = "M = 8", v030 = "settings = \"0.3.0\"")
convergence <- do.call(rbind, lapply(names(arms), function(arm) {
  d <- utils::read.csv(file.path(results, sprintf("panels-K1500-%s.csv", arm)))
  data.frame(
    setting = arms[[arm]], panels = nrow(d), valid = sum(d$valid),
    rhat_max = max(d$max_rhat_parameters), rhat_above_1.01 = sum(d$max_rhat_parameters > 1.01),
    rhat_above_1.05 = sum(d$rhat_above_1.05),
    ess_bulk_min = min(d$min_ess_bulk_parameters), divergences = sum(d$divergences),
    treedepth_share_median = stats::median(d$treedepth_share),
    treedepth_share_max = max(d$treedepth_share),
    leapfrog_median = stats::median(d$leapfrog),
    sampling_seconds_median = stats::median(d$sampling_seconds),
    pilot_seconds_median = stats::median(d$pilot_seconds),
    rmse = mean(d$rmse), inclusion = mean(d$coverage_90),
    inclusion_se = stats::sd(d$coverage_90) / sqrt(nrow(d)),
    mean_length = mean(d$mean_length_90), w1_population = mean(d$w1_population),
    stringsAsFactors = FALSE
  )
}))

# ---- the stored file ---------------------------------------------------------
drop_names <- function(table) {
  rownames(table) <- NULL
  table
}
result <- list(
  krw = list(table = drop_names(krw_table), firms = drop_names(krw_firms),
             cdf = drop_names(krw_cdf)),
  twin = list(table = drop_names(twin_table), cdf = drop_names(twin_cdf)),
  equal_se = list(table = drop_names(equal_table), g = equal_g_table,
                  agreement = equal_agreement, repeated = equal_repeated),
  panels = list(table = panel_table, by_panel = drop_names(panel_rows),
                convergence = drop_names(convergence)),
  settings = data.frame(method = names(calls), call = unname(calls), stringsAsFactors = FALSE),
  about = list(
    created = format(Sys.Date()), interval = "equal-tailed 90%",
    ebnm_intervals = sprintf("sample quantiles of %d draws, seed 1", n_sim),
    seeds = c(krw_fit = 20250807L, twin_towers_fit = 20250627L, range_from_truth = 20250627L,
              equal_se_data = 20261101L, equal_se_fit = 20261102L, panel_data = 20261009L,
              panel_fit = 20261109L, equal_se_repeated_data = 20261200L,
              equal_se_repeated_fit = 20261300L),
    versions = c(
      bayesEfron = as.character(utils::packageVersion("bayesEfron")),
      ebnm = as.character(utils::packageVersion("ebnm")),
      ashr = as.character(utils::packageVersion("ashr")),
      deconvolveR = as.character(utils::packageVersion("deconvolveR")),
      cmdstan = krw_fit$metadata$cmdstan_version
    ),
    protocol = "data-raw/ebnm-protocol.md", script = "data-raw/ebnm-comparison.R",
    notes = c(
      plugged_in = paste(
        "The rows of the g of bayesEfron, plugged in, take the posterior mean of g as known.",
        "Their means, standard deviations and log likelihood are those of ebnm_npmle() with",
        "fix_g = TRUE; their intervals are computed from the same g as the package computes",
        "them for a fit, so that they differ from the rows of bayesEfron by the uncertainty",
        "about g alone."
      ),
      npmle = paste(
        "ebnm_npmle() maximizes the likelihood over a grid that spans the range of the",
        "estimates. A distribution on a wider grid can have a slightly larger likelihood, so",
        "its log likelihood is not an upper limit for the other rows."
      ),
      flat = "A flat prior has no distribution function and no log likelihood.",
      deconvolver = paste(
        "ebnm_deconvolver() takes z-scores. For estimates with unequal standard errors it",
        "estimates the distribution of the effects divided by their standard errors, which",
        "is another quantity, so it is compared only where the standard errors are equal."
      ),
      loglik = paste(
        "The log likelihood of a row of bayesEfron is the marginal log likelihood of the",
        "estimates at its posterior mean of g."
      )
    )
  )
)

# Apart from the times and the date, a second run gives the stored file.
destination <- file.path("inst", "examples", "ebnm_comparison.rds")
without_times <- function(value) {
  value$about$created <- NULL
  for (part in c("krw", "twin", "equal_se")) {
    value[[part]]$table$seconds <- NULL
  }
  value$panels$convergence[c("sampling_seconds_median", "pilot_seconds_median")] <- NULL
  value
}
if (file.exists(destination)) {
  same <- all.equal(without_times(readRDS(destination)), without_times(result), tolerance = 1e-10)
  message("compared with the stored file: ", if (isTRUE(same)) "the same" else paste(same, collapse = "; "))
}
saveRDS(result, destination, compress = "xz", version = 3)

# The two tables of the panels are also kept as text in inst/verification/,
# where the documentation reads its numbers about repeated samples.
stamp <- function(table) {
  cbind(table, generated_on = format(Sys.Date()),
        cmdstan_version = krw_fit$metadata$cmdstan_version, stringsAsFactors = FALSE)
}
utils::write.csv(stamp(result$panels$convergence),
                 file.path("inst", "verification", "panels-convergence.csv"), row.names = FALSE)
utils::write.csv(stamp(result$panels$table),
                 file.path("inst", "verification", "panels-methods.csv"), row.names = FALSE)

if (!is.null(audit_dir)) {
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(table, name) {
    utils::write.csv(table, file.path(audit_dir, name), row.names = FALSE)
  }
  write(result$krw$table, "ebnm-krw.csv")
  write(result$krw$firms, "ebnm-krw-firms.csv")
  write(result$twin$table, "ebnm-twin.csv")
  write(result$equal_se$table, "ebnm-equal-se.csv")
  write(result$equal_se$agreement, "ebnm-equal-se-agreement.csv")
  write(result$equal_se$repeated, "ebnm-equal-se-repeated.csv")
  write(result$panels$table, "ebnm-panels.csv")
  write(result$panels$by_panel, "ebnm-panels-by-panel.csv")
  write(result$panels$convergence, "panels-convergence.csv")
}
options(width = 220, digits = 4)
print(result$krw$table, row.names = FALSE)
print(result$twin$table, row.names = FALSE)
print(result$equal_se$table, row.names = FALSE)
print(result$equal_se$agreement, row.names = FALSE)
print(result$equal_se$repeated, row.names = FALSE)
print(result$panels$table, row.names = FALSE)
print(result$panels$convergence, row.names = FALSE)
message(sprintf("%s: %d bytes; inst/examples/: %.2f MB", destination, file.size(destination),
                sum(file.size(list.files(file.path("inst", "examples"), full.names = TRUE))) / 2^20))
