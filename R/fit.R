#' Fit the Bayesian Efron log-spline model to estimates and standard errors
#'
#' @description
#' A meta-analysis or a multisite trial gives an estimate
#' \eqn{\hat\theta_i}{theta_hat_i} and a standard error
#' \eqn{\sigma_i}{sigma_i} for each of \eqn{K} sites. `bayes_efron_fit()`
#' treats the true effects \eqn{\theta_i}{theta_i} as draws from an
#' unknown distribution \eqn{g} and estimates \eqn{g} and each
#' \eqn{\theta_i}{theta_i} together. Unlike the usual random-effects model,
#' it does not assume that \eqn{g} is normal: \eqn{g} can be skewed or have
#' more than one peak, and the amount of shrinkage for each estimate
#' depends on that shape.
#'
#' The result contains posterior draws of \eqn{g}. Functions such as
#' [site_summary()], [g_summary()] and `plot()` compute posterior means,
#' intervals and probabilities for the sites and for \eqn{g} from these
#' draws.
#'
#' The model is fitted by Markov chain Monte Carlo with CmdStan, so the
#' cmdstanr package and CmdStan must be installed. The first call compiles
#' the Stan program, which can take up to a minute. Later calls reuse the
#' compiled model (see [bayes_efron_compile()]).
#'
#' @details
#' The model for sites \eqn{i = 1, \ldots, K} is
#' \deqn{\hat\theta_i \mid \theta_i \sim N(\theta_i, \sigma_i^2),
#'   \qquad \theta_i \sim g.}{
#'   theta_hat_i | theta_i ~ N(theta_i, sigma_i^2),   theta_i ~ g.}
#' The standard errors \eqn{\sigma_i}{sigma_i} are treated as known and
#' may differ between sites. The distribution \eqn{g} has probabilities
#' \eqn{g_1, \ldots, g_L} on a grid of points
#' \eqn{\tau_1 < \cdots < \tau_L}{tau_1 < ... < tau_L}, with
#' \deqn{\log g_j = B_j^\top \alpha
#'   - \log \sum_{l=1}^{L} \exp(B_l^\top \alpha),}{
#'   log g_j = B_j' alpha - log sum_{l = 1, ..., L} exp(B_l' alpha),}
#' where \eqn{B_j} is the vector of values of the \eqn{M} natural cubic
#' spline functions at \eqn{\tau_j}{tau_j}. This is the log-spline family
#' of Efron (2016). The coefficients have independent priors
#' \eqn{\alpha_m \mid \lambda \sim N(0, 1/\lambda)}{
#' alpha_m | lambda ~ N(0, 1/lambda)},
#' which pull \eqn{g} toward the uniform distribution on the grid. The
#' precision has the prior
#' \eqn{\lambda \sim \mbox{half-Cauchy}(0, 5)}{lambda ~ half-Cauchy(0, 5)},
#' so the strength of that pull is estimated along with the coefficients.
#' See Lee and Sui (2025) for the model and for a comparison with the
#' empirical Bayes estimator, in which \eqn{\alpha}{alpha} is fixed at an
#' estimate. With `component = "normal"`, each grid point carries a normal
#' distribution instead of a point mass, all with the same standard
#' deviation, and \eqn{g} is their mixture with weights \eqn{g_j}.
#'
#' The grid and the spline basis are built as in [make_efron_grid()], which
#' describes the arguments that control them: `grid_method` or `support`
#' for the range, `L` for the number of points, `component` and `width` for
#' what is placed at each point, and `M` and `knots` for the spline
#' functions. Alternatively, a grid built in advance can be supplied as
#' `grid`, in which case the other grid arguments must be omitted.
#'
#' With the default `grid_method = "pilot"`, the range is chosen in two
#' steps. A short preliminary fit, the pilot fit, is run on the range of
#' `"paper_realdata"`, that is, the range of the estimates extended on
#' each side by `expansion` times its width. With the default `expansion`
#' this is the range used in the code that accompanies Lee and Sui (2025).
#' The range is then narrowed to the region where the pilot fit places
#' \eqn{g}. It is never made wider than the initial range, and for some
#' data, e.g., [raudenbush1985] and [krw_discrimination], it is not
#' narrowed at all. See `vignette("choosing-a-grid")` for a discussion of
#' these choices.
#'
#' Only \eqn{\alpha}{alpha} and \eqn{\lambda}{lambda} are sampled. Given
#' \eqn{\alpha}{alpha}, the posterior distribution of each
#' \eqn{\theta_i}{theta_i} has a closed form, and functions such as
#' [site_summary()] compute it from the draws of \eqn{g}. The fitted object
#' stores those draws as one variable per grid point, so its size grows
#' with `L`. With the default of 12,000 draws it takes about 12 MB of
#' memory for 101 grid points and about 40 MB for 401.
#'
#' To make a fit reproducible, specify `seed`. With the same seed, data and
#' settings, the draws are the same on the same machine with the same
#' versions of the packages and of CmdStan. The seed and the sampler
#' settings are stored in the result (see [bef_fit]).
#'
#' A warning is issued if there are divergent transitions, or if an R-hat
#' value exceeds 1.05 (for fits with at least 400 draws in total). A
#' warning is also issued if a quantity is constant within each chain but
#' differs between chains, so that its R-hat is infinite. Use [diagnose()]
#' to see all the diagnostics.
#'
#' @inheritParams make_efron_grid
#' @param theta_hat Numeric vector of effect estimates, one per site, all
#'   on the same scale (e.g., standardized mean differences or log odds
#'   ratios). Names, if present, are used as site labels. Can also be a
#'   data frame or a `bef_data` object that contains the estimates and
#'   their standard errors (see [as_bef_data()]), in which case `sigma` is
#'   omitted. At least 5 sites are required.
#' @param sigma Numeric vector of standard errors of the estimates, the
#'   same length as `theta_hat`. All must be positive.
#' @param ... These dots are for future extensions and must be empty.
#' @param grid A `bef_grid` object created by [make_efron_grid()]. If
#'   `NULL` (the default), the grid is built from the arguments that
#'   follow. If a grid is supplied, `grid_method`, `support`, `L`,
#'   `expansion`, `margin`, `M`, `knots`, `component`, `width`,
#'   `theta_true` and `bound_expansion` must not be specified.
#' @param grid_method Method for choosing the limits of the grid:
#'   `"pilot"` (the default), `"paper_realdata"`, `"paper_simulation"` or
#'   `"paper_sensitivity"`. See the Details of [make_efron_grid()]. With
#'   `"pilot"`, the pilot fit uses this call's `seed` and
#'   `parallel_chains`, so the grid is the same as the one returned by
#'   `make_efron_grid()` with that seed.
#' @param L Number of grid points, an integer from 51 to 401. If `NULL`
#'   (the default), the number is chosen from the data: 101, or more if
#'   101 points are too few for the standard errors or for the spread of
#'   the effects. See the Details of [make_efron_grid()].
#' @param knots Method for placing the knots of the spline functions.
#'   Currently only `"grid"` is available; see the Details of
#'   [make_efron_grid()].
#' @param settings Set to `"0.3.0"` to fit the model as version 0.3.0 of
#'   the package did. The fit then uses that version's Stan program and its
#'   defaults for `grid_method` (`"paper_realdata"`), `L` (101) and
#'   `site_draws` (`TRUE`); arguments that you set yourself are kept. With
#'   `NULL` (the default), the same seed gives different draws than it did
#'   in version 0.3.0, and a message about the changed defaults is shown
#'   once per session. See `vignette("implementation")` for when the draws
#'   of version 0.3.0 are reproduced exactly.
#' @param model_family Model family. Currently only `"RE"`, the model
#'   described in Details, is available.
#' @param chains Number of Markov chains, an integer from 1 to 16.
#' @param parallel_chains Number of chains to run in parallel, an integer
#'   from 1 to `chains`. By default all chains run in parallel.
#' @param iter_warmup Number of warmup iterations per chain.
#' @param iter_sampling Number of sampling (post-warmup) iterations per
#'   chain, at least 1.
#' @param adapt_delta Target acceptance rate of the sampler, a number
#'   between 0 and 1. A higher value leads to a smaller step size, which
#'   can remove divergent transitions but makes sampling slower. With about
#'   20 sites, divergent transitions are common with the default of 0.9,
#'   and 0.99 is a better starting value (see
#'   `vignette("teacher-expectancy")`).
#' @param max_treedepth Maximum tree depth of the sampler, an integer from
#'   1 to 20. Increase it if [diagnose()] reports many iterations that
#'   reached the maximum.
#' @param seed Integer seed for the sampler. If `NULL` (the default), a
#'   seed is taken from the clock. In either case the seed that was used
#'   is stored in the result.
#' @param keep_cmdstan_fit If `TRUE`, the `cmdstanr::CmdStanMCMC` object is
#'   stored in the result as `cmdstan_fit`. It refers to temporary files,
#'   so it is only useful in the session that created it.
#' @param store_grid_quantities If `TRUE`, the draws also include the grid
#'   probabilities `g` and their unnormalized logarithms `log_w`, which
#'   adds `2 * L` variables. This is not needed for the methods, which
#'   compute both from the other draws.
#' @param site_draws If `TRUE`, the fit also stores five variables per
#'   site (`theta_mean`, `theta_sd`, `theta_map`, `theta_rep` and
#'   `log_lik`), which `confint(method = "draws")` requires. If `FALSE`
#'   (the default, except with `settings = "0.3.0"`), only the draws of the
#'   parameters and of \eqn{g} are stored, and all site summaries are
#'   computed from them when they are requested. Both choices give the
#'   same draws of \eqn{g} for the same seed.
#' @param thin Thinning interval, an integer from 1 to `iter_sampling`.
#'   From each chain, one iteration out of every `thin` is stored, starting
#'   with the first. Thinning reduces the size of the result. It does not
#'   affect
#'   the sampler diagnostics reported by [diagnose()] (e.g., R-hat and the
#'   effective sample sizes), which are computed from all iterations. The
#'   summaries and the grid diagnostics are computed from the stored
#'   draws.
#' @param refresh Number of iterations between progress updates from
#'   CmdStan; `0` turns them off. If `NULL` (the default), 200 is used in
#'   an interactive session and 0 otherwise.
#' @param init Initial values for the chains. A number `x` (default 0.5)
#'   draws the initial values uniformly between `-x` and `x` on the
#'   unconstrained scale, i.e., for the spline coefficients and for the
#'   logarithm of the precision. The default is narrower than CmdStan's
#'   default of 2. Can also be a list with one list of values of `alpha`
#'   and `lambda` per chain, or a function that returns such a list of
#'   values.
#' @param show_messages If `TRUE`, the output of CmdStan and the messages
#'   of cmdstanr about the sampler are shown. If `NULL` (the default), they
#'   are shown in an interactive session only. The package's own warnings
#'   and [diagnose()] do not depend on this argument.
#'
#' @return An object of class `c("bef_fit_re", "bef_fit")`, which is a list
#'   with the following components.
#'   \describe{
#'     \item{`draws`}{The posterior draws, a `posterior::draws_array`.}
#'     \item{`posterior`}{A list with the draws of the mean, the variance
#'       and the standard deviation of \eqn{g}, of `effective_params` and
#'       of `log_marginal_likelihood`.}
#'     \item{`metadata`}{A list with the data, the settings and summaries
#'       of the draws.}
#'     \item{`cmdstan_fit`}{The `cmdstanr::CmdStanMCMC` object, only with
#'       `keep_cmdstan_fit = TRUE`.}
#'   }
#'   See [bef_fit] for the details of these components and for the
#'   available methods.
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
#' # Estimates and standard errors as two vectors. The call needs CmdStan.
#' \dontrun{
#' fit <- bayes_efron_fit(
#'   theta_hat = krw_discrimination$yi,
#'   sigma = krw_discrimination$sei,
#'   seed = 20250807
#' )
#'
#' # The call that made `krw_fit`: the firms' names as site labels, and
#' # the draws thinned to every 15th.
#' dat <- as_bef_data(krw_discrimination, site = "firm_name")
#' fit <- bayes_efron_fit(dat, seed = 20250807, thin = 15)
#' }
#'
#' # The stored result of the second call.
#' fit <- krw_fit
#' fit
#' summary(fit)$prior_summary
#' confint(fit, parm = c(29, 65))
#' plot(fit, type = "g")
#'
#' @export
bayes_efron_fit <- function(theta_hat,
                            sigma = NULL,
                            ...,
                            grid = NULL,
                            grid_method = "pilot",
                            support = NULL,
                            L = NULL,
                            expansion = 0.5,
                            margin = 0.1,
                            M = 6L,
                            knots = "grid",
                            component = c("point", "normal"),
                            width = 2 / 3,
                            theta_true = NULL,
                            bound_expansion = NULL,
                            settings = NULL,
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
                            site_draws = FALSE,
                            thin = 1L,
                            refresh = NULL,
                            init = 0.5,
                            show_messages = NULL) {
  call <- match.call()
  settings <- .bef_match_settings(settings)
  given <- setdiff(names(call), "")
  rule_given <- "grid_method" %in% given
  if (!is.null(grid)) {
    .bef_check_grid_alone(intersect(.bef_fit_grid_args(), given))
  } else if (identical(settings, "0.3.0")) {
    # The defaults of version 0.3.0 take the place of those that were not
    # given.
    if (!rule_given) grid_method <- "paper_realdata"
    if (!"L" %in% given) L <- 101L
  }
  if (identical(settings, "0.3.0") && !"site_draws" %in% given) {
    site_draws <- TRUE
  }
  if (is.null(settings) && is.null(grid) && is.null(support) && !rule_given) {
    .bef_inform_default_grid()
  }

  context <- .bef_fit_prepare_context(
    call = call,
    theta_hat = theta_hat,
    sigma = sigma,
    ...,
    grid = grid,
    grid_method = grid_method,
    rule_given = rule_given,
    support = support,
    L = L,
    expansion = expansion,
    margin = margin,
    M = M,
    knots = knots,
    component = component,
    width = width,
    theta_true = theta_true,
    bound_expansion = bound_expansion,
    settings = settings,
    model_family = model_family,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit,
    store_grid_quantities = store_grid_quantities,
    site_draws = site_draws,
    thin = thin,
    refresh = refresh,
    init = init,
    show_messages = show_messages
  )

  context <- .bef_fit_run_sampler(context)
  .bef_fit_assemble(context)
}

# The arguments of bayes_efron_fit() that describe the grid. With `grid`
# none of them may be given.
.bef_fit_grid_args <- function() {
  c("grid_method", "support", "L", "expansion", "margin", "M", "knots", "component",
    "width", "theta_true", "bound_expansion")
}

.bef_check_grid_alone <- function(also) {
  if (length(also) == 0L) {
    return(invisible(TRUE))
  }
  .bef_abort_invalid_args(
    sprintf(
      "`grid` already holds the settings of the grid, so %s must be left out. To change them, build another grid with make_efron_grid().",
      paste(sprintf("`%s`", also), collapse = ", ")
    ),
    arg = "grid"
  )
}

.bef_match_settings <- function(settings) {
  if (is.null(settings)) {
    return(NULL)
  }
  if (!is.character(settings) || length(settings) != 1L || is.na(settings) ||
      !identical(settings, "0.3.0")) {
    .bef_abort_invalid_args(
      "`settings` must be `NULL`, for the defaults of this version, or \"0.3.0\".",
      arg = "settings"
    )
  }
  settings
}

# Said once in a session, at the first fit that relies on the default rule
# for the range of the grid.
.bef_inform_default_grid <- function() {
  key <- "default_grid_0.4"
  if (isTRUE(.bayesEfron_msgs_emitted[[key]])) {
    return(invisible(FALSE))
  }
  .bayesEfron_msgs_emitted[[key]] <- TRUE
  rlang::inform(
    paste(
      "The defaults of bayes_efron_fit() changed in version 0.4: the range of the grid",
      "is set by a short pilot fit (`grid_method = \"pilot\"`) and can be narrower than",
      "before, and the fit keeps the draws of g only. `settings = \"0.3.0\"` gives the",
      "fit of version 0.3.0. This message is shown once in a session."
    ),
    class = "bef_default_changed"
  )
  invisible(TRUE)
}

# The default rule of this function is the one that needs no fit; it is
# what tests of the pieces use. bayes_efron_fit() passes its own.
.bef_fit_prepare_context <- function(call,
                                     theta_hat,
                                     sigma = NULL,
                                     ...,
                                     grid = NULL,
                                     grid_method = "paper_realdata",
                                     rule_given = FALSE,
                                     support = NULL,
                                     L = NULL,
                                     expansion = 0.5,
                                     margin = 0.1,
                                     M = 6L,
                                     knots = "grid",
                                     component = "point",
                                     width = 2 / 3,
                                     theta_true = NULL,
                                     bound_expansion = NULL,
                                     settings = NULL,
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
                                     site_draws = FALSE,
                                     thin = 1L,
                                     refresh = NULL,
                                     init = 0.5,
                                     show_messages = NULL,
                                     check_installed = TRUE,
                                     check_installed_fun = .bef_check_cmdstanr_installed,
                                     now = Sys.time) {
  # The program is chosen by `settings` alone.
  settings <- .bef_match_settings(settings)
  program <- if (identical(settings, "0.3.0")) "efron_re" else "efron_mix"
  given <- .bef_resolve_data_args(theta_hat, sigma)
  args <- .bef_validate_fit_args(
    theta_hat = given$theta_hat,
    sigma = given$sigma,
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
    store_grid_quantities = store_grid_quantities,
    site_draws = site_draws,
    thin = thin,
    refresh = refresh,
    init = init,
    show_messages = show_messages
  )

  if (isTRUE(check_installed)) {
    check_installed_fun()
  }

  # The pilot fit of the rule "pilot" runs with the seed of the fit.
  effective_seed <- .bef_effective_seed(args$seed, now = now)
  grid <- if (!is.null(grid)) {
    .bef_fit_given_grid(grid, theta_hat = args$theta_hat, sigma = args$sigma)
  } else {
    .bef_build_grid(
      theta_hat = args$theta_hat,
      sigma = args$sigma,
      L = args$L,
      expansion = args$expansion,
      M = args$M,
      grid_method = args$grid_method,
      theta_true = args$theta_true,
      bound_expansion = args$bound_expansion,
      support = support,
      component = component,
      width = width,
      knots = knots,
      margin = margin,
      seed = effective_seed,
      rule_given = rule_given,
      parallel_chains = args$parallel_chains
    )
  }
  # The settings as they were used, whoever chose them.
  args$grid_method <- grid$grid_method
  args$L <- grid$L
  args$M <- grid$M
  args$settings <- settings
  bef_data <- as_bef_data(list(theta_hat = args$theta_hat, sigma = args$sigma))
  stan_data <- prepare_stan_data(
    bef_data = bef_data,
    grid = grid,
    model_family = args$model_family,
    store_grid_quantities = isTRUE(args$store_grid_quantities),
    program = program
  )

  structure(
    list(
      call = call,
      args = args,
      effective_seed = effective_seed,
      grid = grid,
      bef_data = bef_data,
      stan_data = stan_data,
      program = program
    ),
    class = "bef_fit_context"
  )
}

# A grid that the caller built: checked, and compared again with the data of
# this fit, which need not be those that it was built from.
.bef_fit_given_grid <- function(grid, theta_hat, sigma) {
  if (!inherits(grid, "bef_grid")) {
    .bef_abort_invalid_args(
      "`grid` must be a grid made by make_efron_grid(), or `NULL`.",
      arg = "grid"
    )
  }
  .bef_validate_grid_return(grid)
  grid$resolution <- .bef_grid_resolution(
    grid$grid, grid$spacing, theta_hat, sigma, component = grid$component,
    width = grid$width, L_chosen = isTRUE(grid$resolution$L_chosen)
  )
  grid$data <- list(theta_hat = as.numeric(theta_hat), sigma = as.numeric(sigma))
  .bef_warn_if_coarse(grid)
  grid
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

# The Stan program of a fit context. A context without one runs
# efron_re.stan.
.bef_context_program <- function(context) {
  if (is.null(context$program)) "efron_re" else context$program
}

.bef_fit_get_model <- function(context, model_fun = .bef_model) {
  # A replacement of the model function that takes the family only is called
  # without the program.
  takes_program <- any(c("program", "...") %in% names(formals(model_fun)))
  tryCatch(
    do.call(model_fun, c(
      list(model_family = context$args$model_family),
      if (takes_program) list(program = .bef_context_program(context))
    )),
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

.bef_fit_refresh <- function(interactive_fun = interactive, refresh = NULL) {
  if (!is.null(refresh)) {
    return(as.integer(refresh))
  }
  if (isTRUE(interactive_fun())) {
    return(200L)
  }
  0L
}

.bef_fit_sampler_settings <- function(context, interactive_fun = interactive) {
  settings <- list(
    chains = context$args$chains,
    parallel_chains = context$args$parallel_chains,
    iter_warmup = context$args$iter_warmup,
    iter_sampling = context$args$iter_sampling,
    adapt_delta = context$args$adapt_delta,
    max_treedepth = context$args$max_treedepth,
    seed = context$effective_seed,
    refresh = .bef_fit_refresh(interactive_fun, context$args$refresh),
    init = if (is.null(context$args$init)) 0.5 else context$args$init
  )
  # Outside an interactive session CmdStan's output and the messages of
  # cmdstanr about the sampler are off unless they are asked for. Thinning
  # is done after sampling, so CmdStan is not told about it.
  show <- context$args$show_messages
  if (is.null(show)) {
    show <- isTRUE(interactive_fun())
  }
  if (!isTRUE(show)) {
    settings$show_messages <- FALSE
    settings$show_exceptions <- FALSE
    settings["diagnostics"] <- list(NULL)
  }
  if (identical(.bef_context_program(context), "efron_mix")) {
    # Every summary of a site is computed from the draws of log_g as CmdStan
    # writes them, and the number of digits that CmdStan writes by default
    # has changed between versions. The program writes few columns, so
    # twelve digits cost nothing. The call for efron_re.stan is left as it
    # was in version 0.3.0.
    settings$sig_figs <- 12L
  }
  settings
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
  # A replacement of the postprocessing function is given `site_draws` only
  # if it takes it, and `program` and the seed only if it takes `program`.
  takes <- names(formals(postprocess_fun))
  keeps <- intersect(c("site_draws", "..."), takes)
  programs <- intersect(c("program", "..."), takes)
  processed <- tryCatch(
    do.call(postprocess_fun, c(
      list(
        cmdstan_fit = context$cmdstan_fit,
        stan_data = context$stan_data,
        model_family = context$args$model_family
      ),
      if (length(keeps) > 0L) list(site_draws = isTRUE(context$args$site_draws)),
      if (length(programs) > 0L) {
        list(program = .bef_context_program(context), seed = context$effective_seed)
      },
      if (any(c("thin", "...") %in% takes) && !is.null(context$args$thin) &&
          context$args$thin > 1L) {
        list(thin = context$args$thin)
      }
    )),
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

  # The labels of the sites travel with the table of site effects.
  labels <- context$bef_data$names
  if (!is.null(labels) && is.data.frame(processed$theta_summary) &&
      length(labels) == nrow(processed$theta_summary)) {
    processed$theta_summary$label <- as.character(labels)
  }

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
  fit <- validate_bef_fit_re(fit)
  .bef_fit_check_grid(fit, grid = context$grid)
}

# Look at the grid of a new fit: store what is found with the fit, and warn
# if a site leans on an end of the grid that no fit on a wider range set.
.bef_fit_check_grid <- function(fit, grid = NULL) {
  found <- tryCatch(
    .bef_grid_diagnostics(
      .bef_engine(fit),
      grid = grid,
      grid_method = fit$metadata$grid_method,
      labels = fit$metadata$theta_summary$label,
      B = fit$metadata$data_list$B
    ),
    error = function(err) NULL
  )
  if (is.null(found)) {
    return(fit)
  }
  attr(fit$metadata, "grid_diagnostics") <- found
  if (any(found$ends$status == "warning")) {
    attr(fit$metadata, "sampler_diagnostics_warned") <- union(
      attr(fit$metadata, "sampler_diagnostics_warned", exact = TRUE), "grid_ends"
    )
    .bef_warn_grid_ends(found)
  }
  fit
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
    stan_file_sha256 = if ("program" %in% names(formals(stan_sha_fun))) {
      stan_sha_fun(context$args$model_family, program = .bef_context_program(context))
    } else {
      stan_sha_fun(context$args$model_family)
    },
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
  # What the fit keeps and when it was made.
  dims <- dim(processed$draws)
  attr(metadata, "storage") <- list(
    site_draws = isTRUE(context$args$site_draws),
    iterations = if (length(dims) == 3L) dims[[1L]] else NA_integer_,
    chains = if (length(dims) == 3L) dims[[2L]] else NA_integer_,
    thin = if (is.null(context$args$thin)) 1L else as.integer(context$args$thin),
    created = format(Sys.Date(), "%Y-%m-%d")
  )
  attr(metadata, "program") <- .bef_context_program(context)
  attr(metadata, "settings") <- context$args$settings
  # The grid as make_efron_grid() returns it: the rule, the support, the
  # components, the knots and the resolution. The estimates that it was
  # built from are left out, since the fit holds its own.
  grid <- context$grid
  grid$data <- NULL
  attr(metadata, "grid") <- grid
  attr(metadata, "nuts") <- processed$nuts
  attr(metadata, "sampler_settings") <- context$sampler_settings
  attr(metadata, "package_version") <- package_version_fun("bayesEfron")
  attr(metadata, "cmdstanr_version") <- package_version_fun("cmdstanr")
  metadata
}

.bef_fit_stan_sha256 <- function(model_family, program = "efron_re") {
  .bef_stan_file_sha256(.bef_stan_file(model_family, program = program))
}

.bef_runtime_package_version <- function(package) {
  as.character(getNamespaceVersion(package))
}
