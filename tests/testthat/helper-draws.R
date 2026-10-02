# Build a synthetic draws array holding the site-level quantities, optional
# scalar series and, when `grid_points` is positive, a density on the grid.
# The values must be the ones used to build `theta_summary`:
# validate_bef_fit_re() checks that theta_summary$mean equals
# colMeans(theta_mean) and that theta_summary$sd is the posterior standard
# deviation computed from theta_mean and theta_sd in the draws.

bef_fixture_draws <- function(theta_map, theta_mean, theta_sd, theta_rep,
                              scalars = list(), grid_points = 0L) {
  S <- nrow(theta_rep)
  K <- ncol(theta_rep)
  cols <- list()
  nms <- character()

  add <- function(values, name) {
    cols[[length(cols) + 1L]] <<- as.numeric(values)
    nms <<- c(nms, name)
  }

  for (entry in list(
    list(m = theta_map, nm = "theta_map"),
    list(m = theta_mean, nm = "theta_mean"),
    list(m = theta_sd, nm = "theta_sd"),
    list(m = theta_rep, nm = "theta_rep")
  )) {
    for (k in seq_len(K)) {
      add(entry$m[, k], sprintf("%s[%d]", entry$nm, k))
    }
  }
  for (nm in names(scalars)) {
    add(scalars[[nm]], nm)
  }
  if (grid_points > 0L) {
    for (l in seq_len(grid_points)) {
      # a plausible positive density that varies across draws, so the
      # grid-level quantile band in the plot payload is properly ordered
      add((l / grid_points) * (1 + seq_len(S) / (10 * S)), sprintf("g[%d]", l))
    }
  }

  array(
    unlist(cols),
    dim = c(S, 1L, length(cols)),
    dimnames = list(NULL, NULL, nms)
  )
}
