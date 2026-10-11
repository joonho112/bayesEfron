# Builds the three fits in data/ with the defaults of this version:
#
#   krw_fit          bayes_efron_fit(as_bef_data(krw_discrimination, site = "firm_name"),
#                                    seed = 20250807, thin = 15)
#   raudenbush_fit   bayes_efron_fit(raudenbush1985, seed = 1985, thin = 15)
#   twin_towers_fit  bayes_efron_fit(twin_towers, M = 8, seed = 20250627, thin = 15)
#
# and inst/examples/twin_towers_comparison.rds, a table that sets the fit of
# the 1,500 simulated sites beside fits on other ranges, among them the
# range around the true effects. Needs CmdStan. Run from the package root:
#
#   Rscript data-raw/example-fits.R [audit-dir] [objects]
#
# `objects` names the fits to build, separated by commas; without it all
# three are built. The table of the simulated sites is built with
# `twin_towers_fit`. The fit of the callback gaps carries the names of the
# firms as the labels of its sites.
#
# Each call runs four chains of 1,000 warmup and 3,000 sampling iterations
# and keeps one draw in fifteen, 800 in all. The diagnostics that the fits
# store are those of all 12,000 draws. Every fit is made twice, once with
# all its draws, so that the summaries of the 800 draws can be compared
# with those of the 12,000; the audit directory receives that comparison
# (example-fits-audit.csv) and the table of the simulated sites
# (twin-towers-comparison.csv). Without an audit directory a temporary one
# is used. The script stops if a stored fit misses the targets below.

args <- commandArgs(trailingOnly = TRUE)
audit_dir <- if (length(args) >= 1L) args[[1L]] else file.path(tempdir(), "bayesEfron-example-fits")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
objects <- c("krw_fit", "raudenbush_fit", "twin_towers_fit")
if (length(args) >= 2L) {
  objects <- intersect(objects, strsplit(args[[2L]], ",", fixed = TRUE)[[1L]])
  stopifnot(length(objects) >= 1L)
}

quiet_fit <- function(...) {
  caught <- list()
  fit <- withCallingHandlers(
    bayes_efron_fit(...),
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (inherits(m, "bef_default_changed")) invokeRestart("muffleMessage")
    }
  )
  attr(fit, "warnings") <- vapply(caught, function(w) {
    sprintf("[%s] %s", class(w)[[1L]], conditionMessage(w))
  }, character(1L))
  fit
}

# The diagnostics of a fit in one row. The columns that end in `_parameters`
# are taken over `alpha`, `lambda`, `mean_g` and `sd_g`, on which the
# convergence target below is checked.
diagnostics_row <- function(fit) {
  d <- attr(fit$metadata, "diagnostics", exact = TRUE)
  core <- c("lambda", "mean_g", "sd_g", grep("^alpha\\[", names(d$rhat), value = TRUE))
  iterations <- attr(fit$metadata, "sampler_settings", exact = TRUE)
  n_run <- iterations$chains * iterations$iter_sampling
  grid <- attr(fit$metadata, "grid", exact = TRUE)
  data.frame(
    L = fit$metadata$data_list$L, M = fit$metadata$data_list$M,
    lower = grid$support[1L], upper = grid$support[2L],
    rule = fit$metadata$grid_method,
    share_of_wide_range = if (is.null(grid$pilot)) NA_real_ else grid$pilot$share,
    max_rhat = max(d$rhat, na.rm = TRUE), max_rhat_parameters = max(d$rhat[core]),
    min_ess_bulk = min(d$ess_bulk, na.rm = TRUE), min_ess_bulk_parameters = min(d$ess_bulk[core]),
    min_ess_tail = min(d$ess_tail, na.rm = TRUE),
    divergences = sum(d$divergences), treedepth_hits = sum(d$max_treedepth),
    treedepth_share = sum(d$max_treedepth) / n_run, min_ebfmi = min(d$ebfmi),
    pilot_seconds = if (is.null(grid$pilot)) NA_real_ else grid$pilot$seconds,
    sampling_seconds = fit$metadata$runtime_seconds,
    mean_g = mean(fit$posterior$mean_g), sd_g = mean(fit$posterior$sd_g),
    stringsAsFactors = FALSE
  )
}

# ---- the three fits ----------------------------------------------------------
cases <- list(
  krw_fit = list(data = as_bef_data(krw_discrimination, site = "firm_name"),
                 data_name = "as_bef_data(krw_discrimination, site = \"firm_name\")",
                 args = list(seed = 20250807L)),
  raudenbush_fit = list(data = raudenbush1985, data_name = "raudenbush1985",
                        args = list(seed = 1985L)),
  twin_towers_fit = list(data = twin_towers, data_name = "twin_towers",
                         args = list(M = 8L, seed = 20250627L))
)
audit <- list()
full_fits <- list()
for (name in objects) {
  case <- cases[[name]]
  full <- do.call(quiet_fit, c(list(case$data), case$args))
  stored <- do.call(quiet_fit, c(list(case$data), case$args, list(thin = 15L)))
  full_fits[[name]] <- full
  stopifnot(
    identical(dim(stored$draws)[1:2], c(200L, 4L)),
    identical(unclass(stored$draws)[, , "lambda"],
              unclass(posterior::thin_draws(full$draws, 15L))[, , "lambda"]),
    identical(attr(stored$metadata, "diagnostics", exact = TRUE),
              attr(full$metadata, "diagnostics", exact = TRUE)),
    is.null(stored$cmdstan_fit),
    !attr(stored$metadata, "storage", exact = TRUE)$site_draws
  )

  # The summaries of the 800 draws beside those of the 12,000.
  spacing <- diff(stored$metadata$data_list$grid[1:2])
  a <- site_summary(full)
  b <- site_summary(stored)
  limits <- abs(as.matrix(b[, c("lower", "upper")]) - as.matrix(a[, c("lower", "upper")])) / spacing
  conditional <- .bef_site_conditional(.bef_engine(stored))$mean
  mcse <- apply(conditional, 2L, stats::sd) / sqrt(nrow(conditional))
  positive <- abs(exceedance_prob(stored, threshold = 0) - exceedance_prob(full, threshold = 0))

  row <- cbind(
    data.frame(object = name, sites = nobs(stored), seed = stored$metadata$seed,
               call = sprintf("bayes_efron_fit(%s%s, seed = %d, thin = 15)", case$data_name,
                              if (is.null(case$args$M)) "" else sprintf(", M = %d", case$args$M),
                              case$args$seed),
               draws_run = posterior::ndraws(full$draws),
               draws_kept = posterior::ndraws(stored$draws),
               stringsAsFactors = FALSE),
    diagnostics_row(stored),
    data.frame(
      warnings = length(attr(stored, "warnings")),
      limit_shift_max_spacings = max(limits),
      limit_shift_p95_spacings = unname(stats::quantile(limits, 0.95)),
      limit_shift_median_spacings = stats::median(limits),
      limits_beyond_two_spacings = sum(limits > 2),
      mean_shift_max = max(abs(b$mean - a$mean)),
      mean_shift_max_mcse = max(abs(b$mean - a$mean) / mcse),
      positive_probability_shift_max = max(positive),
      sd_g_shift = mean(stored$posterior$sd_g) - mean(full$posterior$sd_g)
    )
  )

  # The targets of a stored fit: convergence of the run, no warning, and
  # limits of the site intervals within two spacings of those of the run.
  stopifnot(
    row$max_rhat_parameters < 1.01, row$min_ess_bulk_parameters >= 400,
    row$divergences == 0, row$warnings == 0,
    row$limit_shift_max_spacings <= 2
  )

  attr(stored, "warnings") <- NULL
  assign(name, stored)
  destination <- file.path("data", paste0(name, ".rda"))
  save(list = name, file = destination, compress = "xz", version = 3)
  restored <- new.env(parent = emptyenv())
  load(destination, envir = restored)
  stopifnot(identical(get(name, envir = restored), stored))
  row$bytes <- file.size(destination)
  row$megabytes_in_memory <- as.numeric(utils::object.size(stored)) / 2^20
  audit[[name]] <- row
  message(sprintf("%s: %d bytes, largest shift of an interval limit %.2f spacings",
                  destination, row$bytes, row$limit_shift_max_spacings))
}
audit <- do.call(rbind, audit)
audit$cmdstan_version <- get(objects[[1L]])$metadata$cmdstan_version
audit$created <- format(Sys.Date())
utils::write.csv(audit, file.path(audit_dir, "example-fits-audit.csv"), row.names = FALSE)
options(width = 220)
print(audit[, c("object", "sites", "L", "M", "lower", "upper", "max_rhat", "min_ess_bulk",
                "divergences", "treedepth_hits", "sampling_seconds",
                "limit_shift_max_spacings", "bytes")], row.names = FALSE)
if (!"twin_towers_fit" %in% objects) {
  quit(save = "no")
}

# ---- the simulated sites on other ranges -------------------------------------
# The Wasserstein distance between g, read with its mass spread evenly over
# the cell of each grid point, and the true effects, over an interval that
# holds every range compared here; the root mean squared error of the
# posterior means; and the share of the true effects inside the 90%
# intervals.
truth <- twin_towers$theta_true
wide <- make_efron_grid(twin_towers, grid_method = "paper_realdata")$support
x <- seq(wide[1L] - 1, wide[2L] + 1, length.out = 40001L)
truth_cdf <- stats::ecdf(truth)(x)
recovery <- function(fit) {
  engine <- .bef_engine(fit)
  g <- colMeans(exp(engine$log_g))
  cdf <- as.numeric(g %*% .bef_point_mass_below(engine, x, "interpolated"))
  table <- site_summary(fit)
  data.frame(
    w1 = sum(abs(cdf - truth_cdf)) * (x[2L] - x[1L]),
    rmse = sqrt(mean((table$mean - truth)^2)),
    coverage_90 = mean(truth >= table$lower & truth <= table$upper)
  )
}
specs <- list(
  list(label = "example: rule \"pilot\", M = 8", fit = full_fits$twin_towers_fit),
  list(label = "defaults: rule \"pilot\", M = 6", args = list()),
  list(label = "range around the true effects, M = 6",
       args = list(grid_method = "paper_simulation", theta_true = truth)),
  list(label = "range around the true effects, M = 8",
       args = list(grid_method = "paper_simulation", theta_true = truth, M = 8L)),
  list(label = "range of version 0.3.0 (\"paper_realdata\"), M = 6",
       args = list(grid_method = "paper_realdata"))
)
rows <- lapply(specs, function(spec) {
  fit <- if (is.null(spec$fit)) {
    do.call(quiet_fit, c(list(twin_towers, seed = 20250627L), spec$args))
  } else {
    spec$fit
  }
  message("fitted: ", spec$label)
  cbind(data.frame(fit = spec$label, uses_true_effects = identical(
    fit$metadata$grid_method, "paper_simulation"), stringsAsFactors = FALSE),
    diagnostics_row(fit), recovery(fit),
    data.frame(warnings = paste(attr(fit, "warnings"), collapse = " | "),
               stringsAsFactors = FALSE))
})
comparison <- do.call(rbind, rows)
reference <- comparison[comparison$fit == "range around the true effects, M = 6", ]
comparison$w1_over_reference <- comparison$w1 / reference$w1
comparison$rmse_over_reference <- comparison$rmse / reference$rmse
utils::write.csv(comparison, file.path(audit_dir, "twin-towers-comparison.csv"),
                 row.names = FALSE)
twin_towers_comparison <- list(
  table = comparison,
  reference = "range around the true effects, M = 6",
  seed = 20250627L,
  draws = 12000L,
  cmdstan_version = twin_towers_fit$metadata$cmdstan_version,
  created = format(Sys.Date())
)
saveRDS(twin_towers_comparison, file.path("inst", "examples", "twin_towers_comparison.rds"),
        version = 3)

print(comparison[, c("fit", "lower", "upper", "max_rhat_parameters", "min_ess_bulk_parameters",
                     "divergences", "treedepth_share", "sampling_seconds", "w1", "rmse",
                     "coverage_90", "w1_over_reference", "rmse_over_reference")],
      row.names = FALSE, digits = 4)
message(sprintf("data/: %.2f MB; inst/examples/: %.2f MB",
                sum(file.size(list.files("data", full.names = TRUE))) / 2^20,
                sum(file.size(list.files(file.path("inst", "examples"), full.names = TRUE))) / 2^20))
