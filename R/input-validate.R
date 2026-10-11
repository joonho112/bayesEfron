.bef_validate_fit_args <- function(theta_hat,
                                   sigma,
                                   ...,
                                   grid_method = "paper_realdata",
                                   L = NULL,
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
                                   site_draws = FALSE,
                                   thin = 1L,
                                   refresh = NULL,
                                   init = 0.5,
                                   show_messages = NULL) {
  .bef_check_no_dots(list(...))

  theta_hat <- .bef_validate_theta_hat(theta_hat)
  sigma <- .bef_validate_sigma(sigma, K = length(theta_hat))
  grid_method <- .bef_validate_grid_method(grid_method)
  model_family <- .bef_validate_fit_model_family(model_family)

  # `NULL` leaves the number of grid points to the rule of the grid.
  if (!is.null(L)) {
    L <- .bef_validate_fit_int(
      L,
      arg = "L",
      lower = 51L,
      upper = 401L,
      predicate = "integer scalar in [51, 401], or `NULL`"
    )
  }
  expansion <- .bef_validate_fit_number(
    expansion,
    arg = "expansion",
    lower = 0,
    upper = 5,
    predicate = "finite number in [0, 5]"
  )
  M <- .bef_validate_fit_int(
    M,
    arg = "M",
    lower = 3L,
    upper = 10L,
    predicate = "integer scalar in [3, 10]"
  )
  theta_true <- .bef_validate_theta_true(theta_true, K = length(theta_hat))
  bound_expansion <- .bef_validate_bound_expansion(bound_expansion)

  chains <- .bef_validate_fit_int(
    chains,
    arg = "chains",
    lower = 1L,
    upper = 16L,
    predicate = "integer scalar in [1, 16]"
  )
  parallel_chains <- .bef_validate_fit_int(
    parallel_chains,
    arg = "parallel_chains",
    lower = 1L,
    upper = chains,
    predicate = "integer scalar between 1 and `chains`"
  )
  iter_warmup <- .bef_validate_fit_int(
    iter_warmup,
    arg = "iter_warmup",
    lower = 0L,
    predicate = "non-negative integer scalar"
  )
  iter_sampling <- .bef_validate_fit_int(
    iter_sampling,
    arg = "iter_sampling",
    lower = 1L,
    predicate = "positive integer scalar"
  )
  adapt_delta <- .bef_validate_adapt_delta(adapt_delta)
  max_treedepth <- .bef_validate_fit_int(
    max_treedepth,
    arg = "max_treedepth",
    lower = 1L,
    upper = 20L,
    predicate = "integer scalar in [1, 20]"
  )
  seed <- .bef_validate_seed(seed)
  keep_cmdstan_fit <- .bef_validate_keep_cmdstan_fit(keep_cmdstan_fit)
  store_grid_quantities <- .bef_validate_store_grid_quantities(
    store_grid_quantities
  )
  .bef_check_fit_arg(
    checkmate::assert_flag(site_draws),
    arg = "site_draws",
    predicate = "single TRUE/FALSE value"
  )
  thin <- .bef_validate_fit_int(
    thin,
    arg = "thin",
    lower = 1L,
    upper = iter_sampling,
    predicate = "integer scalar between 1 and `iter_sampling`"
  )
  if (!is.null(refresh)) {
    refresh <- .bef_validate_fit_int(
      refresh,
      arg = "refresh",
      lower = 0L,
      upper = .Machine$integer.max,
      predicate = "non-negative integer scalar, or `NULL`"
    )
  }
  init <- .bef_validate_init(init, chains = chains)
  if (!is.null(show_messages)) {
    .bef_check_fit_arg(
      checkmate::assert_flag(show_messages),
      arg = "show_messages",
      predicate = "single TRUE/FALSE value, or `NULL`"
    )
  }

  list(
    theta_hat = theta_hat,
    sigma = sigma,
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
}

# The starting values of the chains: one number, the half-width of the
# interval from which CmdStan draws them on the unconstrained scale; a list
# with one list of values for each chain; or a function that returns such a
# list of values.
.bef_validate_init <- function(init, chains) {
  number <- is.numeric(init) && length(init) == 1L && is.finite(init) && init >= 0
  per_chain <- is.list(init) && length(init) == chains && all(vapply(init, is.list, logical(1L)))
  if (!number && !per_chain && !is.function(init)) {
    .bef_abort_invalid_args(
      "`init` must be a single number that is not negative, a list with one list of starting values for each chain, or a function that returns such values.",
      arg = "init"
    )
  }
  init
}

.bef_check_no_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }

  dot_names <- names(dots)
  dot_names <- dot_names[nzchar(dot_names)]
  unsupported <- if (length(dot_names) > 0L) {
    paste(sprintf("`%s`", dot_names), collapse = ", ")
  } else {
    "unnamed arguments"
  }

  .bef_abort_invalid_args(
    sprintf("Unused arguments: %s.", unsupported),
    arg = "..."
  )
}

.bef_validate_theta_hat <- function(theta_hat) {
  problem <- .bef_vector_problem(theta_hat, "theta_hat", min_len = .bef_min_sites())
  if (!is.null(problem)) {
    .bef_abort_invalid_args(problem, arg = "theta_hat")
  }
  .bef_drop_storage_mode(theta_hat)
}

.bef_validate_sigma <- function(sigma, K) {
  problem <- .bef_vector_problem(sigma, "sigma", len = K, positive = TRUE)
  if (!is.null(problem)) {
    .bef_abort_invalid_args(problem, arg = "sigma")
  }
  .bef_drop_storage_mode(sigma)
}

# What is wrong with a vector of estimates or standard errors, as one
# sentence, or NULL. The type, missing values, infinite values, the length
# and the sign are reported separately, the first that applies.
.bef_vector_problem <- function(x, arg, min_len = NULL, len = NULL, positive = FALSE) {
  count <- function(n, what) sprintf("%d %s%s", n, what, if (n == 1L) "" else "s")
  if (is.null(x)) {
    return(sprintf("`%s` is missing; give one value for each site.", arg))
  }
  if (!is.numeric(x)) {
    return(sprintf("`%s` must be a numeric vector, not %s.", arg,
                   paste(class(x), collapse = "/")))
  }
  if (anyNA(x)) {
    return(sprintf("`%s` has %s (the first at position %d).", arg,
                   count(sum(is.na(x)), "missing value"), which(is.na(x))[1L]))
  }
  if (any(is.infinite(x))) {
    return(sprintf("`%s` has %s (the first at position %d).", arg,
                   count(sum(is.infinite(x)), "infinite value"), which(is.infinite(x))[1L]))
  }
  if (!is.null(min_len) && length(x) < min_len) {
    return(sprintf("`%s` must have at least %d values, one for each site; it has %d.",
                   arg, min_len, length(x)))
  }
  if (!is.null(len) && length(x) != len) {
    return(sprintf("`%s` must have one value for each of the %d estimates; it has %d.",
                   arg, len, length(x)))
  }
  if (isTRUE(positive) && any(x < .Machine$double.eps)) {
    bad <- x < .Machine$double.eps
    return(sprintf("`%s` must be positive; it has %s that %s not (the first at position %d).",
                   arg, count(sum(bad), "value"), if (sum(bad) == 1L) "is" else "are",
                   which(bad)[1L]))
  }
  NULL
}

.bef_validate_grid_method <- function(grid_method) {
  methods <- .bef_grid_methods()
  if (identical(grid_method, "kl_target_experimental")) {
    .bef_abort_invalid_args(.bef_kl_target_removed(), arg = "grid_method")
  }

  .bef_check_fit_arg(
    checkmate::assert_choice(grid_method, choices = methods),
    arg = "grid_method",
    predicate = sprintf("one of %s", paste(sprintf('"%s"', methods), collapse = ", "))
  )
  grid_method
}

.bef_validate_fit_model_family <- function(model_family) {
  if (!identical(model_family, "RE")) {
    .bef_abort_invalid_args(
      "`model_family` must be \"RE\", the only model that the package fits.",
      arg = "model_family"
    )
  }
  model_family
}

.bef_validate_fit_int <- function(x,
                                  arg,
                                  lower = -Inf,
                                  upper = Inf,
                                  predicate = "integer scalar") {
  .bef_check_fit_arg(
    checkmate::assert_int(x, lower = lower, upper = upper),
    arg = arg,
    predicate = predicate
  )
  as.integer(x)
}

.bef_validate_fit_number <- function(x,
                                     arg,
                                     lower = -Inf,
                                     upper = Inf,
                                     predicate = "finite number") {
  .bef_check_fit_arg(
    checkmate::assert_number(x, lower = lower, upper = upper, finite = TRUE),
    arg = arg,
    predicate = predicate
  )
  as.numeric(x)
}

.bef_validate_theta_true <- function(theta_true, K) {
  if (is.null(theta_true)) {
    return(NULL)
  }

  .bef_check_fit_arg(
    checkmate::assert_numeric(
      theta_true,
      finite = TRUE,
      any.missing = FALSE,
      len = K
    ),
    arg = "theta_true",
    predicate = "finite numeric vector matching theta_hat length"
  )

  .bef_drop_storage_mode(theta_true)
}

.bef_validate_bound_expansion <- function(bound_expansion) {
  if (is.null(bound_expansion)) {
    return(NULL)
  }

  bound_expansion <- .bef_validate_fit_number(
    bound_expansion,
    arg = "bound_expansion",
    lower = 0,
    upper = 5,
    predicate = "finite number in (0, 5]"
  )

  if (bound_expansion <= 0) {
    .bef_abort_invalid_args(
      "`bound_expansion` must be greater than 0.",
      arg = "bound_expansion"
    )
  }

  bound_expansion
}

.bef_validate_adapt_delta <- function(adapt_delta) {
  adapt_delta <- .bef_validate_fit_number(
    adapt_delta,
    arg = "adapt_delta",
    lower = 0,
    upper = 1,
    predicate = "finite number in (0, 1)"
  )

  if (adapt_delta <= 0 || adapt_delta >= 1) {
    .bef_abort_invalid_args(
      "`adapt_delta` must be greater than 0 and less than 1.",
      arg = "adapt_delta"
    )
  }

  adapt_delta
}

.bef_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  .bef_validate_fit_int(
    seed,
    arg = "seed",
    lower = 0L,
    upper = .Machine$integer.max,
    predicate = "NULL or non-negative integer scalar"
  )
}

.bef_validate_store_grid_quantities <- function(store_grid_quantities) {
  .bef_check_fit_arg(
    checkmate::assert_flag(store_grid_quantities),
    arg = "store_grid_quantities",
    predicate = "single TRUE/FALSE value"
  )
  store_grid_quantities
}

.bef_validate_keep_cmdstan_fit <- function(keep_cmdstan_fit) {
  .bef_check_fit_arg(
    checkmate::assert_flag(keep_cmdstan_fit),
    arg = "keep_cmdstan_fit",
    predicate = "single TRUE/FALSE value"
  )
  keep_cmdstan_fit
}

.bef_check_fit_arg <- function(expr, arg, predicate) {
  tryCatch(
    {
      force(expr)
      invisible(TRUE)
    },
    error = function(err) {
      # One sentence: what the argument must be and what was found. The
      # message of the check is not repeated as a parent error.
      found <- sub("^Assertion on '[^']*' failed: ", "", conditionMessage(err))
      found <- sub("[.]$", "", found)
      found <- paste0(tolower(substr(found, 1L, 1L)), substr(found, 2L, nchar(found)))
      .bef_abort_invalid_args(
        sprintf("`%s` must be %s %s: %s.", arg,
                if (grepl("^[aeiou]", predicate)) "an" else "a", predicate, found),
        arg = arg
      )
    }
  )
}

.bef_drop_storage_mode <- function(x) {
  out <- as.numeric(x)
  names(out) <- names(x)
  out
}

# The estimates and standard errors of a call. The first argument is the
# numeric vector of estimates, with `sigma` their standard errors, or a
# `bef_data` object or a data frame that holds both, with `sigma` left out.
# The labels of the sites, when there are any, are returned as the names of
# the estimates.
.bef_resolve_data_args <- function(theta_hat, sigma = NULL) {
  if (!inherits(theta_hat, "bef_data") && !is.data.frame(theta_hat)) {
    return(list(theta_hat = theta_hat, sigma = sigma))
  }
  if (!is.null(sigma)) {
    .bef_abort_invalid_args(
      paste0(
        "`sigma` must be left out when the first argument is a data frame or a ",
        "`bef_data` object, which holds the standard errors."
      ),
      arg = "sigma"
    )
  }
  data <- as_bef_data(theta_hat)
  estimates <- data$theta_hat
  names(estimates) <- data$names
  list(theta_hat = estimates, sigma = data$sigma)
}
