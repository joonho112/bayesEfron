# The range rule "pilot".
#
# The rule starts from the range of "paper_realdata" and narrows it where
# the distribution of the effects occupies only part of it. It never makes
# the range wider.
#
#   (a) The wide range: the range of the estimates, widened on each side by
#       `expansion` times its width.
#   (b) A short fit of the model on that range (the pilot fit) gives a first
#       estimate of g, the posterior mean of its probabilities.
#   (c) The core: the grid points between which the central 99.9% of that
#       estimate lies.
#   (d) The core is widened on each side by `margin` times its width, and by
#       at least one spacing of the pilot grid.
#   (e) If an estimate lies more than four of its standard errors beyond an
#       end of that range, the end is moved out to that estimate plus or
#       minus two standard errors. A site so far outside cannot be explained
#       by an effect inside the range.
#   (f) Neither end goes beyond the wide range.
#
# Parts (b) to (d) and (f), with the settings of the pilot fit below, were
# chosen among several rules by simulation; those simulations are not part
# of this repository. Part (e) and the lower bound in (d) are safeguards:
# they changed none of the 1,759 ranges that the rule narrowed in them.

# The settings of the pilot fit and of the rule. They are not arguments:
# the rule was tested with these values.
.bef_pilot_settings <- function() {
  list(
    L = 101L,
    M = 6L,
    chains = 4L,
    iter_warmup = 500L,
    iter_sampling = 500L,
    adapt_delta = 0.9,
    max_treedepth = 10L,
    # The share of the pilot estimate of g that the core leaves out, half of
    # it at each end.
    eps = 0.001,
    # An estimate more than `conflict` standard errors outside the range
    # moves the end out, to the estimate plus or minus `cover` standard
    # errors.
    conflict = 4,
    cover = 2
  )
}

# The range of "paper_realdata": the range of the estimates, widened on each
# side by `expansion` times its width.
.bef_realdata_range <- function(theta_hat, expansion) {
  range_theta <- range(theta_hat, na.rm = TRUE)
  width <- diff(range_theta)
  if (!is.finite(width) || width <= 0) {
    .bef_abort_grid(
      "The estimates are all equal, so they give no range for the grid. Give the two ends in `support`."
    )
  }
  pad <- expansion * width
  c(range_theta[1L] - pad, range_theta[2L] + pad)
}

# Parts (c) to (f) of the rule: the range from the pilot estimate `g` of the
# probabilities on `grid`, the grid of the pilot fit. `wide` holds the two
# end points of the wide range; where the rule narrows nothing, the result
# has exactly these numbers.
.bef_pilot_support <- function(grid, g, theta_hat, sigma, margin,
                               wide = c(grid[1L], grid[length(grid)]),
                               settings = .bef_pilot_settings()) {
  L <- length(grid)
  spacing <- (wide[2L] - wide[1L]) / (L - 1L)
  if (length(g) != L || any(!is.finite(g)) || any(g < 0) || sum(g) <= 0) {
    .bef_abort_grid(
      "The pilot fit did not return a distribution on its grid.",
      class = "bef_pilot_failed"
    )
  }
  cumulative <- cumsum(g / sum(g))
  first <- which(cumulative >= settings$eps / 2)[1L]
  last <- which(cumulative >= 1 - settings$eps / 2)[1L]
  if (is.na(last)) {
    last <- L
  }
  core <- c(grid[first], grid[last])
  pad <- max(margin * (core[2L] - core[1L]), spacing)
  range <- c(core[1L] - pad, core[2L] + pad)

  below <- (range[1L] - theta_hat) / sigma > settings$conflict
  above <- (theta_hat - range[2L]) / sigma > settings$conflict
  if (any(below)) {
    range[1L] <- min(range[1L], (theta_hat - settings$cover * sigma)[below])
  }
  if (any(above)) {
    range[2L] <- max(range[2L], (theta_hat + settings$cover * sigma)[above])
  }

  support <- c(max(range[1L], wide[1L]), min(range[2L], wide[2L]))
  list(
    support = support,
    wide = wide,
    core = core,
    protected = which(below | above),
    narrowed = support[1L] > wide[1L] || support[2L] < wide[2L],
    share = (support[2L] - support[1L]) / (wide[2L] - wide[1L])
  )
}

# The pilot fit: a short run of the program of this version with point
# masses on `grid`. Returns the posterior mean of the probabilities and what
# the sampler reported. Nothing is printed, and the warnings of CmdStan
# about the sampler are not raised: a pilot fit on a wide range often
# reaches the maximum tree depth, which does not matter for where the bulk
# of g lies.
.bef_pilot_fit <- function(theta_hat, sigma, grid, seed, parallel_chains = NULL,
                           settings = .bef_pilot_settings(),
                           model_fun = .bef_model,
                           check_installed_fun = .bef_check_cmdstanr_installed) {
  need <- paste(
    "The rule \"pilot\" fits the model briefly to set the range of the grid,",
    "which needs CmdStan. Without it, `grid_method = \"paper_realdata\"`",
    "gives the range that the rule starts from, and `support` sets the end",
    "points directly."
  )
  tryCatch(
    check_installed_fun(),
    error = function(err) .bef_abort_pilot(need, parent = err)
  )
  model <- tryCatch(
    model_fun(program = "efron_mix"),
    bef_error = function(err) .bef_abort_pilot(need, parent = err)
  )
  L <- length(grid)
  data <- list(
    K = length(theta_hat),
    theta_hat = as.numeric(theta_hat),
    sigma = as.numeric(sigma),
    L = L,
    grid = as.numeric(grid),
    M = settings$M,
    B = unname(unclass(.bef_grid_basis(grid, M = settings$M, knots = "grid"))),
    tau_c = 0,
    store_grid_quantities = 0L
  )
  attributes(data$B) <- list(dim = c(L, settings$M))
  if (is.null(parallel_chains)) {
    parallel_chains <- settings$chains
  }
  started <- Sys.time()
  fit <- tryCatch(
    model$sample(
      data = data,
      chains = settings$chains,
      parallel_chains = min(as.integer(parallel_chains), settings$chains),
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      adapt_delta = settings$adapt_delta,
      max_treedepth = settings$max_treedepth,
      seed = as.integer(seed),
      refresh = 0L,
      init = 0.5,
      sig_figs = 12L,
      show_messages = FALSE,
      show_exceptions = FALSE,
      diagnostics = NULL
    ),
    interrupt = function(err) stop(err),
    error = function(err) {
      .bef_abort_pilot("CmdStan failed in the pilot fit of the rule \"pilot\".", parent = err)
    }
  )
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  out <- tryCatch(
    {
      log_g <- fit$draws("log_g", format = "draws_matrix")
      log_g <- matrix(as.numeric(log_g), nrow = nrow(log_g))
      if (nrow(log_g) < settings$chains * settings$iter_sampling) {
        stop("some chains did not finish")
      }
      list(g = colMeans(exp(log_g)), diagnostics = .bef_pilot_diagnostics(fit))
    },
    error = function(err) {
      .bef_abort_pilot(
        "The pilot fit of the rule \"pilot\" did not return the draws of g.",
        parent = err
      )
    }
  )
  # The files of the pilot fit are not needed again.
  try(unlink(c(fit$output_files(), fit$data_file())), silent = TRUE)
  out$seconds <- seconds
  out
}

# What the sampler reported for the pilot fit. These numbers are kept with
# the grid; they raise no warning.
.bef_pilot_diagnostics <- function(fit) {
  tryCatch(
    {
      sampler <- fit$diagnostic_summary(quiet = TRUE)
      summary <- fit$summary(c("alpha", "lambda"), "rhat", "ess_bulk")
      list(
        max_rhat = max(summary$rhat),
        min_ess_bulk = min(summary$ess_bulk),
        divergences = sum(sampler$num_divergent),
        max_treedepth = sum(sampler$num_max_treedepth)
      )
    },
    error = function(err) {
      list(max_rhat = NA_real_, min_ess_bulk = NA_real_, divergences = NA_integer_,
           max_treedepth = NA_integer_)
    }
  )
}

.bef_abort_pilot <- function(message, parent = NULL) {
  rlang::abort(
    message,
    class = c("bef_pilot_failed", "bef_grid_error", "bef_error"),
    parent = parent,
    call = .bef_user_frame()
  )
}

# The grid of the rule "pilot". `pilot_fun` is replaced in tests.
.bef_grid_pilot <- function(theta_hat, sigma, expansion, margin, seed,
                            parallel_chains = NULL, pilot_fun = .bef_pilot_fit) {
  settings <- .bef_pilot_settings()
  wide <- .bef_realdata_range(theta_hat, expansion)
  pilot_grid <- seq(from = wide[1L], to = wide[2L], length.out = settings$L)
  seed <- .bef_effective_seed(seed)
  attribution <- list(
    formula = paste(
      "range(theta_hat) widened on each side by expansion * width, narrowed to the",
      "central 99.9% of a pilot estimate of g widened on each side by margin * width;",
      "never wider than the first range"
    ),
    source = "bayesEfron 0.4.0, rule \"pilot\""
  )
  # A range narrower than the median standard error is not narrowed
  # further, and no fit is run: the data cannot place g within it.
  if (diff(wide) < .bef_resolution_rule()$width * stats::median(sigma)) {
    return(list(
      support = wide, expansion = expansion, attribution = attribution,
      pilot = list(
        wide = wide, core = wide, margin = margin, narrowed = FALSE, share = 1,
        protected = integer(), seed = seed, chains = 0L, iter_warmup = 0L,
        iter_sampling = 0L,
        diagnostics = list(max_rhat = NA_real_, min_ess_bulk = NA_real_,
                           divergences = NA_integer_, max_treedepth = NA_integer_),
        seconds = 0
      )
    ))
  }
  pilot <- pilot_fun(
    theta_hat = theta_hat, sigma = sigma, grid = pilot_grid, seed = seed,
    parallel_chains = parallel_chains
  )
  found <- .bef_pilot_support(
    grid = pilot_grid, g = pilot$g, theta_hat = theta_hat, sigma = sigma,
    margin = margin, wide = wide, settings = settings
  )
  # Where the rule narrows nothing, the end points are those of
  # "paper_realdata" exactly, and so is the grid.
  list(
    support = found$support,
    expansion = expansion,
    attribution = attribution,
    pilot = list(
      wide = wide,
      core = found$core,
      margin = margin,
      narrowed = found$narrowed,
      share = found$share,
      protected = found$protected,
      seed = seed,
      chains = settings$chains,
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      diagnostics = pilot$diagnostics,
      seconds = pilot$seconds
    )
  )
}
