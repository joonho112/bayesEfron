#' Fit the Bayesian Efron random-effects model
#'
#' @description
#' Fit the fully Bayesian Efron log-spline prior to a univariate
#' random-effects meta-analytic deconvolution problem. Given per-site
#' effect estimates and their within-study standard errors,
#' `bayes_efron_fit()` returns posterior site-effect summaries together
#' with a continuous estimate of the underlying mixing distribution.
#'
#' The function targets applied meta-analysts who already have point
#' estimates and standard errors on a comparable scale (for example,
#' the `yi` and `sqrt(vi)` columns of an [metafor::escalc()] object).
#' One call performs the full eight-stage pipeline: input validation,
#' grid construction, Stan-data preparation, cache-backed model
#' retrieval, NUTS sampling, draw extraction, generated-quantity
#' postprocessing, and assembly into a validated `bef_fit_re` object.
#'
#' Computation is delegated to CmdStan through `cmdstanr`; the first
#' call in a session typically spends most of its time compiling the
#' Stan program. Subsequent calls reuse the on-disk and in-session
#' caches as long as the model source, package version, CmdStan
#' version, and platform match.
#'
#' @details
#' # The fitted four-level hierarchy
#'
#' For sites \eqn{i = 1, \ldots, K},
#' \deqn{\hat\theta_i \mid \theta_i \sim \mathcal{N}(\theta_i, \sigma_i^2),}
#' \deqn{\theta_i \mid g \overset{\text{iid}}{\sim} g,}
#' where \eqn{g} is the mixing distribution of latent site effects.
#' The package places a log-spline prior on the discretized \eqn{g}
#' with a half-Cauchy hyperprior on the smoothness precision
#' \eqn{\lambda}; full derivations are in the methodological vignettes.
#'
#' # Grid recipes
#'
#' The four `grid_method` choices control how the discrete support of
#' \eqn{g} is constructed:
#'
#' | Recipe | Needs `theta_true`? | Use when |
#' |:-------|:-------------------:|:---------|
#' | `"paper_realdata"` | No | Real-data analysis with no oracle. |
#' | `"paper_simulation"` | Yes | Simulation with known truth, matched paper. |
#' | `"paper_sensitivity"` | Yes | Bound-expansion sensitivity, paper rule. |
#' | `"kl_target_experimental"` | No | KL-target tuning (experimental, heteroscedastic). |
#'
#' The two oracle-requiring recipes refuse to run without a numeric
#' `theta_true`. The experimental recipe emits a once-per-session
#' disclaimer about its KL calibration.
#'
#' # Sampler defaults
#'
#' v0.1 fixes the CmdStan initialization at `init = 0.5` to keep
#' release fits reproducible; this is **not** a user-facing tuning
#' argument. `parallel_chains` defaults to `chains` but can be
#' overridden by setting the environment variable
#' `BAYESEFRON_PARALLEL_CHAINS` (or the unprefixed `PARALLEL_CHAINS`)
#' to a positive integer that does not exceed `chains`.
#'
#' # Reproducibility
#'
#' If `seed` is `NULL`, the function auto-generates an integer seed
#' from the current time and stores it on `fit$metadata$seed` so the
#' fit can be re-played later. Pass an explicit `seed` for fully
#' reproducible runs.
#'
#' # The 13-field metadata contract
#'
#' `fit$metadata` always has exactly the 13 named fields listed in
#' \strong{Value} below. Four additional payloads are stored as
#' attributes of `fit$metadata` and surfaced through [summary.bef_fit()] and
#' [diagnose()]: `sd_g_summary`, `diagnostics`, `diagnostic_skipped`,
#' and `sampler_diagnostics_failed`.
#'
#' @param theta_hat Numeric vector of per-site effect estimates on a
#'   common scale (e.g. mean differences, log odds ratios).
#' @param sigma Numeric vector of strictly positive per-site standard
#'   errors, on the same scale as `theta_hat` and the same length.
#' @param ... Reserved for future expansion; must be empty in v0.1.
#' @param grid_method Character grid recipe. One of
#'   `"paper_realdata"` (default), `"paper_simulation"`,
#'   `"paper_sensitivity"`, or `"kl_target_experimental"`. See
#'   \strong{Details}.
#' @param L Integer grid length (number of discrete support points).
#'   Defaults to `101L`. The package's verification ledger is
#'   calibrated at this default.
#' @param expansion Numeric, non-negative. Range-relative expansion
#'   factor that widens the grid endpoints beyond the observed range
#'   of `theta_hat`. Defaults to `0.5` (50 percent expansion).
#' @param M Integer natural-cubic-spline degrees of freedom. Defaults
#'   to `6L`. The verification ledger is calibrated at this default.
#' @param theta_true Numeric oracle vector of latent site effects,
#'   required by `"paper_simulation"` and `"paper_sensitivity"` and
#'   ignored by the other recipes. Same length as `theta_hat`.
#' @param bound_expansion Numeric, oracle-bound expansion factor used
#'   only by `"paper_sensitivity"`. `NULL` falls back to the recipe
#'   default of `0.5`.
#' @param model_family Character scalar. v0.1 supports `"RE"` only.
#'   Other model families are deferred to v0.2+ per the package
#'   blueprint.
#' @param chains Integer number of MCMC chains. Defaults to `4L`.
#' @param iter_warmup Integer warmup iterations per chain. Defaults
#'   to `1000L`.
#' @param iter_sampling Integer post-warmup iterations per chain.
#'   Defaults to `3000L`.
#' @param adapt_delta NUTS target acceptance statistic in `(0, 1)`.
#'   Defaults to `0.9`.
#' @param seed Integer seed or `NULL`. If `NULL`, an integer seed is
#'   auto-generated and recorded on `fit$metadata$seed`.
#' @param keep_cmdstan_fit Logical. If `TRUE`, the raw
#'   `cmdstanr::CmdStanMCMC` handle is retained at `fit$cmdstan_fit`
#'   for advanced use. Defaults to `FALSE` so the returned object is
#'   small enough to save and share.
#'
#' @return An S3 object of class `c("bef_fit_re", "bef_fit")` with
#'   three top-level fields:
#'
#'   * `draws` — a `posterior::draws_array` of MCMC draws for the
#'     model parameters and generated quantities.
#'   * `metadata` — a named list with exactly 13 fields:
#'     - `model_family` — `"RE"` for v0.1.
#'     - `grid_method` — the recipe used.
#'     - `seed` — effective integer seed (auto-generated if not
#'        supplied).
#'     - `cmdstan_version` — CmdStan version string.
#'     - `stan_file_sha256` — SHA-256 of the locked Stan source.
#'     - `data_list` — the seven-field Stan data block sent to
#'        CmdStan.
#'     - `runtime_seconds` — sampler wall-clock.
#'     - `mean_g_summary`, `var_g_summary` — posterior summaries of
#'        functionals of the mixing distribution \eqn{g}.
#'     - `theta_summary` — posterior summaries of latent site
#'        effects \eqn{\theta_i}.
#'     - `theta_rep_draws` — replicated effects for posterior
#'        predictive checks.
#'     - `effective_params_summary`, `log_marginal_likelihood_summary`
#'        — model-quality summaries.
#'   * `posterior` — a tidy posterior representation used by the S3
#'     methods.
#'
#'   Four additional durable payloads are stored as attributes of
#'   `fit$metadata`:  `sd_g_summary`, `diagnostics`,
#'   `diagnostic_skipped`, and `sampler_diagnostics_failed`. Access
#'   them through [summary.bef_fit()] and [diagnose()] rather than directly.
#'
#'   When `keep_cmdstan_fit = TRUE`, `fit$cmdstan_fit` carries the
#'   raw `cmdstanr::CmdStanMCMC` handle.
#'
#' @seealso
#'   * [as_bef_data()] for converting [metafor::escalc()] objects and
#'     plain lists into the package input class.
#'   * [make_efron_grid()] for building grids outside the fitting
#'     pipeline (for example, to inspect a recipe before fitting).
#'   * [bayes_efron_compile()] for pre-warming the Stan model cache.
#'   * [summary.bef_fit()], [confint.bef_fit_re()], [diagnose()], [plot.bef_fit_re()]
#'     for inspecting the returned object.
#'
#' @examples
#' \dontrun{
#' # Five-site smoke fit. Requires a working CmdStan installation.
#' theta_hat <- c(-0.21, 0.04, 0.19, 0.38, 0.61)
#' sigma     <- c( 0.18, 0.15, 0.22, 0.19, 0.24)
#'
#' fit <- bayes_efron_fit(
#'   theta_hat     = theta_hat,
#'   sigma         = sigma,
#'   L             = 51L,
#'   M             = 3L,
#'   chains        = 1L,
#'   iter_warmup   = 150L,
#'   iter_sampling = 4L,
#'   seed          = 1234L
#' )
#'
#' summary(fit)
#' confint(fit, type = "theta")
#' diagnose(fit)
#' plot(fit, type = "caterpillar")
#' }
#'
#' @export
bayes_efron_fit <- function(theta_hat,
                            sigma,
                            ...,
                            grid_method = c(
                              "paper_realdata",
                              "paper_simulation",
                              "paper_sensitivity",
                              "kl_target_experimental"
                            ),
                            L = 101L,
                            expansion = 0.5,
                            M = 6L,
                            theta_true = NULL,
                            bound_expansion = NULL,
                            model_family = "RE",
                            chains = 4L,
                            iter_warmup = 1000L,
                            iter_sampling = 3000L,
                            adapt_delta = 0.9,
                            seed = NULL,
                            keep_cmdstan_fit = FALSE) {
  call <- match.call()
  context <- .bef_fit_prepare_context(
    call = call,
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    grid_method = grid_method,
    L = L,
    expansion = expansion,
    M = M,
    theta_true = theta_true,
    bound_expansion = bound_expansion,
    model_family = model_family,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit
  )

  context <- .bef_fit_run_stages_6_7(context)
  .bef_fit_assemble(context)
}

.bef_fit_prepare_context <- function(call,
                                     theta_hat,
                                     sigma,
                                     ...,
                                     grid_method = .bef_grid_methods(),
                                     L = 101L,
                                     expansion = 0.5,
                                     M = 6L,
                                     theta_true = NULL,
                                     bound_expansion = NULL,
                                     model_family = "RE",
                                     chains = 4L,
                                     iter_warmup = 1000L,
                                     iter_sampling = 3000L,
                                     adapt_delta = 0.9,
                                     seed = NULL,
                                     keep_cmdstan_fit = FALSE,
                                     check_installed = TRUE,
                                     check_installed_fun = .bef_check_cmdstanr_installed,
                                     now = Sys.time) {
  args <- .bef_validate_fit_args(
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    grid_method = grid_method,
    L = L,
    expansion = expansion,
    M = M,
    theta_true = theta_true,
    bound_expansion = bound_expansion,
    model_family = model_family,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit
  )

  if (isTRUE(check_installed)) {
    check_installed_fun()
  }

  grid <- make_efron_grid(
    theta_hat = args$theta_hat,
    sigma = args$sigma,
    L = args$L,
    expansion = args$expansion,
    M = args$M,
    grid_method = args$grid_method,
    theta_true = args$theta_true,
    bound_expansion = args$bound_expansion
  )
  bef_data <- as_bef_data(list(theta_hat = args$theta_hat, sigma = args$sigma))
  stan_data <- prepare_stan_data(
    bef_data = bef_data,
    grid = grid,
    model_family = args$model_family
  )

  structure(
    list(
      call = call,
      args = args,
      effective_seed = .bef_effective_seed(args$seed, now = now),
      grid = grid,
      bef_data = bef_data,
      stan_data = stan_data
    ),
    class = "bef_fit_context"
  )
}

.bef_effective_seed <- function(seed, now = Sys.time) {
  if (!is.null(seed)) {
    return(as.integer(seed))
  }

  as.integer(as.numeric(now())) %% .Machine$integer.max
}

.bef_fit_run_stages_6_7 <- function(context,
                                    model_fun = .bef_model,
                                    sample_fun = NULL,
                                    now = Sys.time,
                                    interactive_fun = interactive) {
  model <- .bef_fit_get_model(context, model_fun = model_fun)

  started <- now()
  cmdstan_fit <- .bef_fit_sample(
    model = model,
    context = context,
    sample_fun = sample_fun,
    interactive_fun = interactive_fun
  )
  runtime_seconds <- as.numeric(difftime(now(), started, units = "secs"))

  context$model <- model
  context$cmdstan_fit <- cmdstan_fit
  context$runtime_seconds <- runtime_seconds
  context
}

.bef_fit_get_model <- function(context, model_fun = .bef_model) {
  tryCatch(
    model_fun(
      model_name = context$args$model_family,
      check_installed = FALSE
    ),
    error = function(err) {
      if (inherits(err, "bef_error")) {
        stop(err)
      }
      .bef_abort_compile_failed(
        "Failed to retrieve the cached bayesEfron Stan model.",
        stage = 6L,
        model_family = context$args$model_family,
        parent = err,
        call = context$call
      )
    }
  )
}

.bef_fit_sample <- function(model,
                            context,
                            sample_fun = NULL,
                            interactive_fun = interactive) {
  if (is.null(sample_fun)) {
    sample_fun <- tryCatch(model$sample, error = function(err) NULL)
  }
  if (!is.function(sample_fun)) {
    .bef_abort_sampling_failed(
      "Compiled bayesEfron model does not expose a callable `sample()` method.",
      stage = 7L,
      model_family = context$args$model_family,
      call = context$call
    )
  }

  set.seed(context$effective_seed)
  cmdstan_fit <- tryCatch(
    sample_fun(
      data = context$stan_data,
      chains = context$args$chains,
      parallel_chains = .bef_fit_parallel_chains(context),
      iter_warmup = context$args$iter_warmup,
      iter_sampling = context$args$iter_sampling,
      adapt_delta = context$args$adapt_delta,
      seed = context$effective_seed,
      refresh = .bef_fit_refresh(interactive_fun),
      init = 0.5
    ),
    interrupt = function(err) {
      stop(err)
    },
    error = function(err) {
      .bef_abort_sampling_failed(
        "CmdStan failed while sampling the bayesEfron model.",
        stage = 7L,
        model_family = context$args$model_family,
        parent = err,
        call = context$call
      )
    }
  )

  .bef_validate_sampling_completion(cmdstan_fit, context)
  cmdstan_fit
}

.bef_fit_parallel_chains <- function(context, getenv = Sys.getenv) {
  raw <- getenv("BAYESEFRON_PARALLEL_CHAINS", unset = "")
  env_name <- "BAYESEFRON_PARALLEL_CHAINS"
  if (!nzchar(raw)) {
    raw <- getenv("PARALLEL_CHAINS", unset = "")
    env_name <- "PARALLEL_CHAINS"
  }
  if (!nzchar(raw)) {
    return(context$args$chains)
  }

  parallel_chains <- suppressWarnings(as.integer(raw))
  if (
    length(parallel_chains) != 1L ||
      is.na(parallel_chains) ||
      !identical(as.character(parallel_chains), raw) ||
      parallel_chains < 1L ||
      parallel_chains > context$args$chains
  ) {
    .bef_abort_sampling_failed(
      sprintf(
        "`%s` must be a single integer between 1 and the requested chain count.",
        env_name
      ),
      stage = 7L,
      model_family = context$args$model_family,
      call = context$call
    )
  }

  parallel_chains
}

.bef_fit_refresh <- function(interactive_fun = interactive) {
  if (isTRUE(interactive_fun())) {
    return(200L)
  }
  0L
}

.bef_validate_sampling_completion <- function(cmdstan_fit, context) {
  completed_fun <- tryCatch(
    cmdstan_fit$num_chains_completed,
    error = function(err) NULL
  )
  if (!is.function(completed_fun)) {
    return(invisible(TRUE))
  }

  completed <- tryCatch(
    completed_fun(),
    error = function(err) {
      .bef_abort_sampling_failed(
        "Failed to inspect completed CmdStan chains after sampling.",
        stage = 7L,
        model_family = context$args$model_family,
        parent = err,
        call = context$call
      )
    }
  )
  completed <- suppressWarnings(as.integer(completed))
  if (length(completed) != 1L || is.na(completed)) {
    .bef_abort_sampling_failed(
      "CmdStan returned an invalid completed-chain count after sampling.",
      stage = 7L,
      model_family = context$args$model_family,
      chains_completed = completed,
      call = context$call
    )
  }

  requested <- context$args$chains
  if (completed <= 0L) {
    .bef_abort_sampling_failed(
      "CmdStan sampling completed zero chains.",
      stage = 7L,
      model_family = context$args$model_family,
      chains_requested = requested,
      chains_completed = completed,
      call = context$call
    )
  }
  if (completed < requested) {
    .bef_abort_sampling_partial(
      "CmdStan sampling completed fewer chains than requested.",
      stage = 7L,
      model_family = context$args$model_family,
      chains_requested = requested,
      chains_completed = completed,
      call = context$call
    )
  }

  invisible(TRUE)
}

.bef_fit_assemble <- function(context,
                              postprocess_fun = postprocess_stan_draws,
                              cmdstan_version_fun = .bef_cmdstan_version,
                              stan_sha_fun = .bef_fit_stan_sha256) {
  processed <- tryCatch(
    postprocess_fun(
      cmdstan_fit = context$cmdstan_fit,
      stan_data = context$stan_data,
      model_family = context$args$model_family
    ),
    error = function(err) {
      if (inherits(err, "bef_error")) {
        stop(err)
      }
      .bef_abort_extraction_failed(
        "Failed to postprocess bayesEfron CmdStan draws.",
        stage = 8L,
        model_family = context$args$model_family,
        parent = err,
        call = context$call
      )
    }
  )

  metadata <- .bef_fit_metadata(
    context = context,
    processed = processed,
    cmdstan_version_fun = cmdstan_version_fun,
    stan_sha_fun = stan_sha_fun
  )
  cmdstan_fit <- if (isTRUE(context$args$keep_cmdstan_fit)) {
    context$cmdstan_fit
  } else {
    NULL
  }

  fit <- new_bef_fit_re(
    draws = processed$draws,
    metadata = metadata,
    posterior = processed$posterior,
    cmdstan_fit = cmdstan_fit
  )
  validate_bef_fit_re(fit)
}

.bef_fit_metadata <- function(context,
                              processed,
                              cmdstan_version_fun,
                              stan_sha_fun) {
  metadata <- list(
    model_family = context$args$model_family,
    grid_method = context$args$grid_method,
    seed = context$effective_seed,
    cmdstan_version = as.character(cmdstan_version_fun()),
    stan_file_sha256 = stan_sha_fun(context$args$model_family),
    data_list = context$stan_data,
    runtime_seconds = context$runtime_seconds,
    mean_g_summary = processed$mean_g_summary,
    var_g_summary = processed$var_g_summary,
    theta_summary = processed$theta_summary,
    theta_rep_draws = processed$theta_rep_draws,
    effective_params_summary = processed$effective_params_summary,
    log_marginal_likelihood_summary = processed$log_marginal_likelihood_summary
  )

  attr(metadata, "sd_g_summary") <- processed$sd_g_summary
  attr(metadata, "diagnostics") <- processed$diagnostics
  attr(metadata, "diagnostic_skipped") <- processed$diagnostic_skipped
  attr(metadata, "sampler_diagnostics_failed") <-
    processed$sampler_diagnostics_failed
  metadata
}

.bef_fit_stan_sha256 <- function(model_family) {
  .bef_stan_file_sha256(.bef_stan_file(model_family))
}

#' Pre-compile the bayesEfron Stan model
#'
#' @description
#' Pre-warm the CmdStan compilation cache used by
#' [bayes_efron_fit()]. The function compiles (or reattaches a cached
#' build of) the locked v0.1 random-effects Stan model, stores the
#' resulting `cmdstanr::CmdStanModel` in the in-session cache, and
#' runs a tiny post-compile smoke check to confirm the binary is
#' callable.
#'
#' Calling `bayes_efron_compile()` before the first
#' `bayes_efron_fit()` shifts the (potentially long) compilation cost
#' out of the fit pipeline, which is useful when fitting
#' interactively or under a wall-clock budget. It is otherwise
#' optional: `bayes_efron_fit()` triggers the same cache mechanism on
#' first use.
#'
#' @details
#' The cache lives at the location given by the environment variable
#' `BAYESEFRON_CACHE_ROOT` (with a sensible per-user default if not
#' set). The lookup key combines the Stan source SHA-256, the package
#' version, the CmdStan version, and the platform, so a cached binary
#' is only reused when all of those match. The post-compile smoke
#' check uses the same fixed internal sampler initialization as
#' [bayes_efron_fit()] (`init = 0.5`).
#'
#' This function requires the `cmdstanr` package and a working
#' CmdStan toolchain. Both are listed in `Suggests:` rather than
#' `Imports:` so the package can be installed and documented without
#' them.
#'
#' @param model_family Character scalar. v0.1 supports `"RE"` only.
#' @param quiet Logical. If `TRUE` (default), compile and smoke-check
#'   output is suppressed during cache warming.
#' @param force_recompile Logical. If `TRUE`, bypasses cached
#'   artifacts and recompiles from Stan source. Defaults to `FALSE`.
#' @param seed_for_check Non-negative integer scalar used for the
#'   synthetic post-compile smoke check. Defaults to `42L`.
#'
#' @return Invisibly, a `cmdstanr::CmdStanModel` reference attached to
#'   the cache entry. The return value is rarely used directly; it is
#'   returned to support advanced workflows that want to drive the
#'   model object outside the package's pipeline.
#'
#' @seealso
#'   * [bayes_efron_fit()] for the user-facing fit pipeline.
#'   * [bayes_efron_clear_cache()] for cache maintenance and stale
#'     lock recovery.
#'
#' @examples
#' \dontrun{
#' # Pre-warm the cache so the next fit skips the compile cost.
#' bayes_efron_compile()
#' }
#'
#' @export
bayes_efron_compile <- function(model_family = "RE",
                                quiet = TRUE,
                                force_recompile = FALSE,
                                seed_for_check = 42L) {
  .bef_check_cmdstanr_installed()
  .bef_compile_entry(
    model_family = model_family,
    quiet = quiet,
    force_recompile = force_recompile,
    seed_for_check = seed_for_check,
    check_installed = FALSE
  )
}

.bef_compile_entry <- function(model_family = "RE",
                               quiet = TRUE,
                               force_recompile = FALSE,
                               seed_for_check = 42L,
                               model_fun = .bef_model,
                               smoke_check_fun = .bef_compile_smoke_check,
                               check_installed = TRUE) {
  if (isTRUE(check_installed)) {
    .bef_check_cmdstanr_installed()
  }

  model_family <- .bef_validate_compile_model_family(model_family)
  quiet <- .bef_validate_compile_flag(quiet, "quiet")
  force_recompile <- .bef_validate_compile_flag(
    force_recompile, "force_recompile"
  )
  seed_for_check <- .bef_validate_compile_seed(seed_for_check)

  cmdstan_model_fun <- .bef_compile_cmdstan_model_fun(quiet)
  model <- .bef_compile_quietly(
    model_fun(
      model_name = model_family,
      force_recompile = force_recompile,
      cmdstan_model_fun = cmdstan_model_fun,
      check_installed = FALSE
    ),
    quiet = quiet
  )

  .bef_compile_quietly(
    smoke_check_fun(
      model = model,
      model_family = model_family,
      seed_for_check = seed_for_check,
      quiet = quiet
    ),
    quiet = quiet
  )

  invisible(model)
}

.bef_compile_cmdstan_model_fun <- function(quiet) {
  force(quiet)
  function(...) {
    .bef_cmdstan_model(..., quiet = quiet)
  }
}

.bef_compile_quietly <- function(expr, quiet) {
  if (!isTRUE(quiet)) {
    return(force(expr))
  }

  value <- NULL
  utils::capture.output(
    value <- withCallingHandlers(
      force(expr),
      message = function(msg) {
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )
  value
}

.bef_compile_smoke_check <- function(model,
                                     model_family,
                                     seed_for_check,
                                     quiet) {
  sample_fun <- tryCatch(model$sample, error = function(err) NULL)
  if (!is.function(sample_fun)) {
    .bef_abort_compile_failed(
      "Compiled bayesEfron model does not expose a callable `sample()` method.",
      model_family = model_family,
      seed_for_check = seed_for_check
    )
  }

  stan_data <- .bef_compile_smoke_stan_data(model_family)
  output_dir <- tempfile("bayesefron-compile-check-")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  check_fit <- NULL
  tryCatch(
    check_fit <- sample_fun(
      data = stan_data,
      chains = 1L,
      parallel_chains = 1L,
      iter_warmup = 2L,
      iter_sampling = 2L,
      seed = seed_for_check,
      refresh = if (isTRUE(quiet)) 0L else 1L,
      init = 0.5,
      output_dir = output_dir,
      show_messages = !isTRUE(quiet),
      diagnostics = c("divergences", "treedepth")
    ),
    error = function(err) {
      .bef_abort_compile_failed(
        "Post-compile smoke check failed for the bayesEfron Stan model.",
        model_family = model_family,
        seed_for_check = seed_for_check,
        parent = err
      )
    }
  )
  .bef_validate_compile_smoke_fit(
    check_fit,
    model_family = model_family,
    seed_for_check = seed_for_check
  )

  invisible(model)
}

.bef_validate_compile_smoke_fit <- function(check_fit,
                                            model_family,
                                            seed_for_check) {
  completed_fun <- tryCatch(
    check_fit$num_chains_completed,
    error = function(err) NULL
  )
  if (!is.function(completed_fun)) {
    return(invisible(TRUE))
  }

  completed <- tryCatch(completed_fun(), error = function(err) NA_integer_)
  if (!identical(as.integer(completed), 1L)) {
    .bef_abort_compile_failed(
      "Post-compile smoke check did not complete its single chain.",
      model_family = model_family,
      seed_for_check = seed_for_check,
      chains_completed = completed
    )
  }

  invisible(TRUE)
}

.bef_compile_smoke_stan_data <- function(model_family = "RE") {
  theta_hat <- c(-0.45, -0.1, 0, 0.25, 0.55)
  sigma <- c(0.12, 0.18, 0.15, 0.22, 0.2)
  bef_data <- as_bef_data(list(theta_hat = theta_hat, sigma = sigma))
  grid <- make_efron_grid(
    theta_hat = theta_hat,
    sigma = sigma,
    L = 51L,
    expansion = 0.5,
    M = 3L,
    grid_method = "paper_realdata"
  )
  prepare_stan_data(bef_data, grid, model_family = model_family)
}

.bef_validate_compile_model_family <- function(model_family) {
  .bef_check_compile_arg(
    checkmate::assert_choice(model_family, choices = "RE"),
    arg = "model_family",
    predicate = '"RE"'
  )
  model_family
}

.bef_validate_compile_flag <- function(x, arg) {
  .bef_check_compile_arg(
    checkmate::assert_flag(x),
    arg = arg,
    predicate = "single TRUE/FALSE value"
  )
  isTRUE(x)
}

.bef_validate_compile_seed <- function(seed_for_check) {
  .bef_check_compile_arg(
    checkmate::assert_int(
      seed_for_check,
      lower = 0L,
      upper = .Machine$integer.max
    ),
    arg = "seed_for_check",
    predicate = "non-negative integer scalar"
  )
  as.integer(seed_for_check)
}

.bef_check_compile_arg <- function(expr, arg, predicate) {
  tryCatch(
    {
      force(expr)
      invisible(TRUE)
    },
    error = function(err) {
      .bef_abort_invalid_args(
        sprintf(
          "`%s` failed validation (%s): %s",
          arg,
          predicate,
          conditionMessage(err)
        ),
        arg = arg,
        predicate = predicate,
        module = "compile-entry",
        stage = 6L,
        parent = err
      )
    }
  )
}
