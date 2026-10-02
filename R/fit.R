#' Fit the Efron log-spline model to estimates and standard errors
#'
#' @description
#' A meta-analysis or a multisite trial gives an estimate
#' \eqn{\hat\theta_i} and a standard error \eqn{\sigma_i} for each of
#' \eqn{K} sites. `bayes_efron_fit()` treats the true effects
#' \eqn{\theta_i} as draws from an unknown distribution \eqn{g} and
#' estimates \eqn{g} and each \eqn{\theta_i} together. Unlike the usual
#' random-effects model, it does not take \eqn{g} to be normal: \eqn{g}
#' can be skewed or have more than one mode, and how far each estimate is
#' shrunk depends on that shape.
#'
#' The function returns posterior draws of \eqn{g}. For each site and each
#' draw of \eqn{g} it also returns the mean, standard deviation and mode of
#' \eqn{\theta_i} given that \eqn{g}, and one draw of \eqn{\theta_i}. The
#' model is fitted by Markov chain Monte Carlo with CmdStan, so the cmdstanr
#' package and CmdStan must be installed. The first call compiles the Stan
#' model, which can take up to a minute; later calls reuse the compiled
#' model (see [bayes_efron_compile()]).
#'
#' @details
#' For sites \eqn{i = 1, \ldots, K},
#' \deqn{\hat\theta_i \mid \theta_i \sim N(\theta_i, \sigma_i^2),
#'   \qquad \theta_i \sim g.}
#' The standard errors \eqn{\sigma_i} are taken as known and may differ
#' between sites. The distribution \eqn{g} has probabilities
#' \eqn{g_1, \ldots, g_L} on a grid of points
#' \eqn{\tau_1 < \cdots < \tau_L}, with
#' \deqn{\log g_j = B_j^\top \alpha - \log \sum_{l=1}^{L} \exp(B_l^\top \alpha),}
#' where \eqn{B_j} holds the values of the \eqn{M} natural cubic spline
#' functions at \eqn{\tau_j}. This is the log-spline family of Efron
#' (2016). The coefficients have independent priors
#' \eqn{\alpha_m \mid \lambda \sim N(0, 1/\lambda)}, which draw \eqn{g}
#' toward the uniform distribution on the grid, and the precision has the
#' prior \eqn{\lambda \sim \mbox{half-Cauchy}(0, 5)}, so that the strength
#' of that pull is estimated along with the coefficients. Lee and Sui
#' (2025) describe the model and compare it with the empirical Bayes
#' estimator, in which \eqn{\alpha} is fixed at an estimate.
#'
#' The grid and the spline basis are those of [make_efron_grid()], and
#' `grid_method`, `L`, `expansion` and `M` control them. The defaults, 101
#' grid points and six spline functions on a grid that extends beyond the
#' estimates by half their range on each side, follow the application in
#' Lee and Sui (2025, Section 7.2) and its replication code.
#' `vignette("choosing-a-grid")` discusses when to change them.
#'
#' Only \eqn{\alpha} and \eqn{\lambda} are sampled. Given \eqn{\alpha}, the
#' posterior distribution of each \eqn{\theta_i} on the grid has a closed
#' form, and its mean, standard deviation and mode and one draw from it are
#' computed for every posterior draw of \eqn{\alpha}.
#'
#' Each chain starts from spline coefficients drawn uniformly between -0.5
#' and 0.5 and from a precision whose logarithm is drawn in the same way,
#' a narrower interval than CmdStan's default of -2 to 2. Pass a `seed` to
#' make a fit reproducible: with the same data, settings, versions of the
#' packages and of CmdStan, and machine, the draws are the same. The seed
#' and the settings of the sampler are stored with the result; [bef_fit]
#' says where.
#'
#' The function gives a warning if there are divergent transitions and, for
#' a fit with at least 400 draws in all, if an R-hat value exceeds 1.05;
#' [diagnose()] returns all the diagnostics.
#'
#' CmdStan writes the draws with a limited number of significant digits,
#' eight in recent versions. That is ample for effect sizes. Estimates that
#' are very large compared with their standard errors, by a factor of
#' 100,000 or more, should be centered before they are passed, because the
#' saved draws would otherwise lose the digits in which the sites differ.
#'
#' @inheritParams make_efron_grid
#' @param theta_hat Numeric vector of effect estimates, one for each site,
#'   all on the same scale (standardized mean differences or log odds
#'   ratios, for example). At least five are needed.
#' @param sigma Numeric vector of the standard errors of the estimates,
#'   positive and of the same length as `theta_hat`.
#' @param ... These dots are for future extensions and must be empty.
#' @param L Whole number between 51 and 300, the number of grid points.
#'   Defaults to `101L`.
#' @param model_family Character string. `"RE"`, the model described in
#'   Details, is the only choice.
#' @param chains Whole number between 1 and 16, the number of Markov
#'   chains. Defaults to `4L`.
#' @param parallel_chains Whole number between 1 and `chains`, the number
#'   of chains to run at the same time. Defaults to `chains`.
#' @param iter_warmup Whole number, the number of warmup iterations in each
#'   chain. Defaults to `1000L`.
#' @param iter_sampling Whole number of at least 1, the number of iterations
#'   in each chain that are kept. Defaults to `3000L`.
#' @param adapt_delta Number between 0 and 1, the target acceptance rate of
#'   the sampler. Defaults to `0.9`. A value nearer 1 gives smaller steps,
#'   which can remove divergent transitions at the cost of a longer run.
#' @param max_treedepth Whole number between 1 and 20, the maximum tree
#'   depth of the sampler. Defaults to `10L`. Raise it if [diagnose()]
#'   reports many iterations that reached the maximum.
#' @param seed Whole number, the seed of the sampler, or `NULL` (the
#'   default), in which case a seed is taken from the clock. The seed that
#'   was used is stored in the result.
#' @param keep_cmdstan_fit Logical. If `TRUE`, the `cmdstanr::CmdStanMCMC`
#'   object is kept as the element `cmdstan_fit` of the result. It refers to
#'   temporary files, so it is of use only in the session that made it.
#'   Defaults to `FALSE`.
#' @param store_grid_quantities Logical. If `TRUE`, the draws also hold the
#'   grid probabilities `g` and their unnormalized logarithms `log_w`,
#'   which adds `2 * L` variables. Defaults to `FALSE`; the methods compute
#'   both from the other draws when they need them.
#'
#' @return An object of class `bef_fit_re`, which inherits from `bef_fit`.
#'   Its elements are `draws`, the posterior draws as a
#'   `posterior::draws_array`; `posterior`, the draws of the mean, variance
#'   and standard deviation of \eqn{g} and of two summaries of the fit; and
#'   `metadata`, the data, the settings and summaries of the draws.
#'   [bef_fit] describes the elements and the methods.
#'
#' @references
#' Efron, B. (2016). Empirical Bayes deconvolution estimates. *Biometrika*,
#' 103(1), 1-20. \doi{10.1093/biomet/asv068}
#'
#' Lee, J. and Sui, D. (2025). Fully Bayesian inference for meta-analytic
#' deconvolution using Efron's log-spline prior. *Mathematics*, 13(16),
#' 2639. \doi{10.3390/math13162639}
#'
#' @seealso [as_bef_data()], [make_efron_grid()], [bef_fit], [diagnose()],
#'   [plot.bef_fit_re()], `vignette("bayesEfron")`
#'
#' @examples
#' # The call that produced `raudenbush_fit`. It needs CmdStan.
#' \dontrun{
#' fit <- bayes_efron_fit(
#'   theta_hat = raudenbush1985$yi,
#'   sigma = raudenbush1985$sei,
#'   seed = 1985
#' )
#' }
#'
#' # The stored result of that call.
#' fit <- raudenbush_fit
#' fit
#' summary(fit)
#' confint(fit, parm = c(4, 10))
#' plot(fit, type = "g")
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
                            parallel_chains = chains,
                            iter_warmup = 1000L,
                            iter_sampling = 3000L,
                            adapt_delta = 0.9,
                            max_treedepth = 10L,
                            seed = NULL,
                            keep_cmdstan_fit = FALSE,
                            store_grid_quantities = FALSE) {
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
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit,
    store_grid_quantities = store_grid_quantities
  )

  context <- .bef_fit_run_sampler(context)
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
                                     parallel_chains = chains,
                                     iter_warmup = 1000L,
                                     iter_sampling = 3000L,
                                     adapt_delta = 0.9,
                                     max_treedepth = 10L,
                                     seed = NULL,
                                     keep_cmdstan_fit = FALSE,
                                     store_grid_quantities = FALSE,
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
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit,
    store_grid_quantities = store_grid_quantities
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
    model_family = args$model_family,
    store_grid_quantities = isTRUE(args$store_grid_quantities)
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

.bef_fit_run_sampler <- function(context,
                                    model_fun = .bef_model,
                                    sample_fun = NULL,
                                    now = Sys.time,
                                    interactive_fun = interactive) {
  model <- .bef_fit_get_model(context, model_fun = model_fun)
  context$sampler_settings <- .bef_fit_sampler_settings(context, interactive_fun)

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
    model_fun(model_family = context$args$model_family),
    error = function(err) {
      if (inherits(err, "bef_error")) {
        stop(err)
      }
      .bef_abort_compile_failed(
        "Failed to compile or load the bayesEfron Stan model.",
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
      model_family = context$args$model_family,
      call = context$call
    )
  }

  settings <- context$sampler_settings
  if (is.null(settings)) {
    settings <- .bef_fit_sampler_settings(context, interactive_fun)
  }
  cmdstan_fit <- tryCatch(
    do.call(sample_fun, c(list(data = context$stan_data), settings)),
    interrupt = function(err) {
      stop(err)
    },
    error = function(err) {
      .bef_abort_sampling_failed(
        "CmdStan failed while sampling the bayesEfron model.",
        model_family = context$args$model_family,
        parent = err,
        call = context$call
      )
    }
  )

  .bef_validate_sampling_completion(cmdstan_fit, context)
  cmdstan_fit
}

.bef_fit_refresh <- function(interactive_fun = interactive) {
  if (isTRUE(interactive_fun())) {
    return(200L)
  }
  0L
}

.bef_fit_sampler_settings <- function(context, interactive_fun = interactive) {
  list(
    chains = context$args$chains,
    parallel_chains = context$args$parallel_chains,
    iter_warmup = context$args$iter_warmup,
    iter_sampling = context$args$iter_sampling,
    adapt_delta = context$args$adapt_delta,
    max_treedepth = context$args$max_treedepth,
    seed = context$effective_seed,
    refresh = .bef_fit_refresh(interactive_fun),
    init = 0.5
  )
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
      model_family = context$args$model_family,
      chains_completed = completed,
      call = context$call
    )
  }

  requested <- context$args$chains
  if (completed <= 0L) {
    .bef_abort_sampling_failed(
      "CmdStan sampling completed zero chains.",
      model_family = context$args$model_family,
      chains_requested = requested,
      chains_completed = completed,
      call = context$call
    )
  }
  if (completed < requested) {
    .bef_abort_sampling_partial(
      "CmdStan sampling completed fewer chains than requested.",
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
                              stan_sha_fun = .bef_fit_stan_sha256,
                              package_version_fun = .bef_runtime_package_version) {
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
    stan_sha_fun = stan_sha_fun,
    package_version_fun = package_version_fun
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
                              stan_sha_fun,
                              package_version_fun = .bef_runtime_package_version) {
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
    effective_params_summary = processed$effective_params_summary,
    log_marginal_likelihood_summary = processed$log_marginal_likelihood_summary
  )

  attr(metadata, "sd_g_summary") <- processed$sd_g_summary
  attr(metadata, "diagnostics") <- processed$diagnostics
  attr(metadata, "diagnostic_skipped") <- processed$diagnostic_skipped
  attr(metadata, "sampler_diagnostics_failed") <-
    processed$sampler_diagnostics_failed
  attr(metadata, "sampler_diagnostics_warned") <-
    processed$sampler_diagnostics_warned
  attr(metadata, "sampler_settings") <- context$sampler_settings
  attr(metadata, "package_version") <- package_version_fun("bayesEfron")
  attr(metadata, "cmdstanr_version") <- package_version_fun("cmdstanr")
  metadata
}

.bef_fit_stan_sha256 <- function(model_family) {
  .bef_stan_file_sha256(.bef_stan_file(model_family))
}

.bef_runtime_package_version <- function(package) {
  as.character(getNamespaceVersion(package))
}
