#' Prepare Stan Data for the Efron Random-Effects Model
#'
#' `prepare_stan_data()` maps validated package-side data and a grid object to
#' the exact seven-field data block consumed by `inst/stan/efron_re.stan`.
#'
#' @param bef_data A list-like object with numeric `theta_hat` and `sigma`
#'   fields.
#' @param grid A grid list returned by [make_efron_grid()].
#' @param model_family Character scalar. v0.1 supports `"RE"` only.
#' @param group Reserved for future families; must be `NULL` in v0.1.
#' @param rho Reserved for future families; must be `NULL` in v0.1.
#'
#' @return A plain named list with fields `K`, `theta_hat`, `sigma`, `L`,
#'   `grid`, `M`, and `B`.
#' @keywords internal
#' @noRd
prepare_stan_data <- function(bef_data,
                              grid,
                              model_family = "RE",
                              group = NULL,
                              rho = NULL) {
  model_family <- .bef_validate_model_family(model_family)
  if (!identical(model_family, "RE")) {
    .bef_abort_data_prep("`model_family` must be \"RE\" for bayesEfron v0.1.")
  }
  if (!is.null(group)) {
    .bef_abort_data_prep("`group` is reserved for future model families and must be NULL in v0.1.")
  }
  if (!is.null(rho)) {
    .bef_abort_data_prep("`rho` is reserved for future model families and must be NULL in v0.1.")
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

  list(
    K = as.integer(length(theta_hat)),
    theta_hat = as.numeric(theta_hat),
    sigma = as.numeric(sigma),
    L = as.integer(grid_obj$L),
    grid = as.numeric(grid_obj$grid),
    M = as.integer(grid_obj$M),
    B = matrix(as.numeric(grid_obj$B), nrow = grid_obj$L, ncol = grid_obj$M)
  )
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
