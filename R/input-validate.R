.bef_validate_fit_args <- function(theta_hat,
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
                                   keep_cmdstan_fit = FALSE) {
  .bef_check_no_dots(list(...))

  theta_hat <- .bef_validate_theta_hat(theta_hat)
  sigma <- .bef_validate_sigma(sigma, K = length(theta_hat))
  grid_method <- .bef_validate_grid_method(grid_method)
  model_family <- .bef_validate_fit_model_family(model_family)

  L <- .bef_validate_fit_int(
    L,
    arg = "L",
    lower = 51L,
    upper = 300L,
    predicate = "integer scalar in [51, 300]"
  )
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
  iter_warmup <- .bef_validate_fit_int(
    iter_warmup,
    arg = "iter_warmup",
    lower = 0L,
    predicate = "non-negative integer scalar"
  )
  iter_sampling <- .bef_validate_fit_int(
    iter_sampling,
    arg = "iter_sampling",
    lower = 0L,
    predicate = "non-negative integer scalar"
  )
  adapt_delta <- .bef_validate_adapt_delta(adapt_delta)
  seed <- .bef_validate_seed(seed)
  keep_cmdstan_fit <- .bef_validate_keep_cmdstan_fit(keep_cmdstan_fit)

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
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    seed = seed,
    keep_cmdstan_fit = keep_cmdstan_fit
  )
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
    sprintf("`...` is closed in bayesEfron v0.1; unsupported arguments: %s.", unsupported),
    arg = "...",
    predicate = "empty dots",
    module = "input-validate",
    stage = 2L
  )
}

.bef_validate_theta_hat <- function(theta_hat) {
  .bef_check_fit_arg(
    checkmate::assert_numeric(
      theta_hat,
      finite = TRUE,
      any.missing = FALSE,
      min.len = 2L
    ),
    arg = "theta_hat",
    predicate = "finite numeric vector length >= 2"
  )

  .bef_drop_storage_mode(theta_hat)
}

.bef_validate_sigma <- function(sigma, K) {
  .bef_check_fit_arg(
    checkmate::assert_numeric(
      sigma,
      finite = TRUE,
      any.missing = FALSE,
      len = K,
      lower = .Machine$double.eps
    ),
    arg = "sigma",
    predicate = "finite positive numeric vector matching theta_hat length"
  )

  .bef_drop_storage_mode(sigma)
}

.bef_validate_grid_method <- function(grid_method) {
  methods <- .bef_grid_methods()
  if (identical(grid_method, methods)) {
    return(methods[[1L]])
  }

  .bef_check_fit_arg(
    checkmate::assert_choice(grid_method, choices = methods),
    arg = "grid_method",
    predicate = sprintf("one of %s", paste(sprintf('"%s"', methods), collapse = ", "))
  )
  grid_method
}

.bef_validate_fit_model_family <- function(model_family) {
  .bef_check_fit_arg(
    checkmate::assert_choice(model_family, choices = "RE"),
    arg = "model_family",
    predicate = '"RE"'
  )
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
      arg = "bound_expansion",
      predicate = "finite number in (0, 5]",
      module = "input-validate",
      stage = 2L
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
      arg = "adapt_delta",
      predicate = "finite number in (0, 1)",
      module = "input-validate",
      stage = 2L
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
      .bef_abort_invalid_args(
        sprintf(
          "`%s` failed validation (%s): %s",
          arg,
          predicate,
          conditionMessage(err)
        ),
        arg = arg,
        predicate = predicate,
        module = "input-validate",
        stage = 2L,
        parent = err
      )
    }
  )
}

.bef_drop_storage_mode <- function(x) {
  out <- as.numeric(x)
  names(out) <- names(x)
  out
}
