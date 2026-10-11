# The data that a Stan program of the package reads, from the estimates and
# standard errors in `bef_data` and a grid from make_efron_grid(): a list
# with K, theta_hat, sigma, L, grid, M, B and store_grid_quantities, and for
# the program "efron_mix" also tau_c, the standard deviation of the component
# at each grid point.
prepare_stan_data <- function(bef_data,
                              grid,
                              model_family = "RE",
                              store_grid_quantities = FALSE,
                              program = "efron_re") {
  model_family <- .bef_validate_model_family(model_family)
  if (!identical(model_family, "RE")) {
    .bef_abort_data_prep("`model_family` must be \"RE\"; other families are not implemented.")
  }

  theta_hat <- .bef_extract_numeric_field(bef_data, "theta_hat", min_len = 2L)
  sigma <- .bef_extract_numeric_field(bef_data, "sigma", min_len = 2L)
  if (length(theta_hat) != length(sigma)) {
    .bef_abort_data_prep("`theta_hat` and `sigma` must have the same length.")
  }
  if (any(sigma <= 0)) {
    .bef_abort_data_prep("`sigma` must be strictly positive before Stan data preparation.")
  }

  grid_obj <- .bef_validate_grid_for_stan(grid)

  if (!is.logical(store_grid_quantities) ||
      length(store_grid_quantities) != 1L ||
      is.na(store_grid_quantities)) {
    .bef_abort_data_prep("`store_grid_quantities` must be a single non-missing logical.")
  }

  program <- match.arg(program, .bef_stan_programs())
  tau_c <- if (is.null(grid$tau_c)) 0 else grid$tau_c
  if (!is.numeric(tau_c) || length(tau_c) != 1L || !is.finite(tau_c) || tau_c < 0) {
    .bef_abort_data_prep("`grid$tau_c` must be a single number that is not negative.")
  }
  if (identical(program, "efron_re") && tau_c > 0) {
    .bef_abort_data_prep(
      "The Stan program of version 0.3.0 has point masses only; it cannot fit a grid with normal components."
    )
  }

  data <- list(
    K = as.integer(length(theta_hat)),
    theta_hat = as.numeric(theta_hat),
    sigma = as.numeric(sigma),
    L = as.integer(grid_obj$L),
    grid = as.numeric(grid_obj$grid),
    M = as.integer(grid_obj$M),
    B = matrix(as.numeric(grid_obj$B), nrow = grid_obj$L, ncol = grid_obj$M)
  )
  if (identical(program, "efron_mix")) {
    data$tau_c <- as.numeric(tau_c)
  }
  data$store_grid_quantities <- as.integer(store_grid_quantities)
  data
}

.bef_validate_model_family <- function(model_family) {
  if (!is.character(model_family) || length(model_family) != 1L ||
      is.na(model_family)) {
    .bef_abort_data_prep("`model_family` must be a single character value.")
  }
  model_family
}

.bef_extract_numeric_field <- function(x, field, min_len) {
  if (!is.list(x) || is.null(x[[field]])) {
    .bef_abort_data_prep(sprintf("`bef_data` must contain `%s`.", field))
  }
  value <- x[[field]]
  if (!is.numeric(value) || length(value) < min_len || any(!is.finite(value))) {
    .bef_abort_data_prep(
      sprintf(
        "`bef_data$%s` must be a finite numeric vector of length at least %d.",
        field, min_len
      )
    )
  }
  as.numeric(value)
}

.bef_validate_grid_for_stan <- function(grid) {
  required <- c("grid", "B", "M", "L")
  if (!is.list(grid) || !all(required %in% names(grid))) {
    .bef_abort_data_prep(
      "`grid` must be a list with fields `grid`, `B`, `M`, and `L`."
    )
  }

  support <- grid$grid
  L <- grid$L
  M <- grid$M
  B <- grid$B

  if (!.bef_is_integerish_scalar(L) || L < 1L) {
    .bef_abort_data_prep("`grid$L` must be a positive integer scalar.")
  }
  if (!.bef_is_integerish_scalar(M) || M < 1L) {
    .bef_abort_data_prep("`grid$M` must be a positive integer scalar.")
  }
  L <- as.integer(L)
  M <- as.integer(M)

  if (!is.numeric(support) || length(support) != L ||
      any(!is.finite(support)) || any(diff(support) <= 0)) {
    .bef_abort_data_prep(
      "`grid$grid` must be a finite strictly increasing numeric vector of length `grid$L`."
    )
  }
  if (!is.matrix(B) || !is.numeric(B) || nrow(B) != L || ncol(B) != M ||
      any(!is.finite(B))) {
    .bef_abort_data_prep(
      "`grid$B` must be a finite numeric matrix with dimensions `grid$L` by `grid$M`."
    )
  }

  list(
    grid = as.numeric(support),
    B = B,
    M = M,
    L = L
  )
}

.bef_is_integerish_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == as.integer(x)
}

.bef_abort_data_prep <- function(message, class = "bef_invalid_args") {
  rlang::abort(
    message,
    class = c(class, "bef_pipeline_error", "bef_error")
  )
}
