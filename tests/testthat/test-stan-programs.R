# The two Stan programs of the package: the program of version 0.3.0, which
# must stay as it was, and the program of this version.

programs_ns <- function(name) getFromNamespace(name, "bayesEfron")

# The hash of a text file with its line ends written as "\n", so that a
# checkout with other line ends gives the same value.
programs_text_sha256 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  digest::digest(paste0(paste(lines, collapse = "\n"), "\n"), algo = "sha256",
                 serialize = FALSE)
}

test_that("the program of version 0.3.0 is unchanged", {
  ref <- readRDS(test_path("_fixtures", "ref-v030.rds"))
  expect_identical(
    programs_text_sha256(programs_ns(".bef_stan_file")(program = "efron_re")),
    ref$conditions$stan_file_sha256
  )
})

test_that("the package has two programs, and an unknown name is an error", {
  expect_identical(programs_ns(".bef_stan_programs")(), c("efron_mix", "efron_re"))
  for (program in programs_ns(".bef_stan_programs")()) {
    file <- programs_ns(".bef_stan_file")(program = program)
    expect_true(file.exists(file))
    expect_identical(basename(file), paste0(program, ".stan"))
  }
  expect_identical(basename(programs_ns(".bef_stan_file")()), "efron_re.stan")
  expect_error(programs_ns(".bef_stan_file")(program = "efron"))
})

test_that("the new program reads tau_c and writes nothing for the sites", {
  code <- readLines(programs_ns(".bef_stan_file")(program = "efron_mix"))
  expect_true(any(grepl("real<lower=0> tau_c;", code, fixed = TRUE)))
  expect_false(any(grepl("_rng", code, fixed = TRUE)))
  expect_false(any(grepl("theta_(map|mean|sd|rep)|effective_params|vector\\[K\\] log_lik;", code)))
})

test_that("the data of the two programs differ by tau_c only", {
  dat <- as_bef_data(list(theta_hat = c(-0.2, 0, 0.3, 0.5, 0.6),
                          sigma = c(0.1, 0.2, 0.15, 0.1, 0.2)))
  grid <- make_efron_grid(dat$theta_hat, dat$sigma, L = 51L, M = 3L, grid_method = "paper_realdata")
  old <- prepare_stan_data(dat, grid)
  new <- prepare_stan_data(dat, grid, program = "efron_mix")
  expect_named(old, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "store_grid_quantities"))
  expect_named(new, c("K", "theta_hat", "sigma", "L", "grid", "M", "B", "tau_c",
                      "store_grid_quantities"))
  expect_identical(new[names(old)], old)
  expect_identical(new$tau_c, 0)

  normal <- make_efron_grid(dat$theta_hat, dat$sigma, L = 51L, M = 3L,
                            grid_method = "paper_realdata", component = "normal")
  expect_gt(normal$tau_c, 0)
  expect_identical(prepare_stan_data(dat, normal, program = "efron_mix")$tau_c, normal$tau_c)
  err <- tryCatch(prepare_stan_data(dat, normal), error = identity)
  expect_s3_class(err, "bef_error")
  expect_match(conditionMessage(err), "point masses only")

  for (bad in list(-0.1, NA_real_, c(0.1, 0.2), "0.1")) {
    broken <- normal
    broken$tau_c <- bad
    expect_s3_class(
      tryCatch(prepare_stan_data(dat, broken, program = "efron_mix"), error = identity),
      "bef_error"
    )
  }
})

# Draws as a plain matrix: one row per draw, chain after chain, and one named
# column per variable.
programs_plain <- function(draws) {
  draws <- posterior::as_draws_matrix(draws)
  matrix(as.numeric(draws), nrow = nrow(draws), dimnames = list(NULL, colnames(draws)))
}

# A stand-in for the object that CmdStan returns, holding what efron_mix.stan
# writes for given draws of the coefficients.
programs_fake_fit <- function(alpha, data, n_chains = 2L) {
  n_draws <- nrow(alpha)
  n_iter <- n_draws %/% n_chains
  log_w <- alpha %*% t(data$B)
  log_g <- log_w - apply(log_w, 1L, function(x) max(x) + log(sum(exp(x - max(x)))))
  g <- exp(log_g)
  mean_g <- drop(g %*% data$grid)
  var_g <- rowSums(g * outer(mean_g, data$grid, function(m, x) (x - m)^2)) + data$tau_c^2
  s_tot <- sqrt(data$sigma^2 + data$tau_c^2)
  lml <- vapply(seq_len(n_draws), function(s) {
    sum(vapply(seq_len(data$K), function(i) {
      log(sum(g[s, ] * stats::dnorm(data$theta_hat[[i]], data$grid, s_tot[[i]])))
    }, numeric(1L)))
  }, numeric(1L))
  values <- cbind(
    lp__ = -seq_len(n_draws) / 10, alpha, lambda = seq(0.5, 1.5, length.out = n_draws),
    log_g, mean_g = mean_g, var_g = var_g, sd_g = sqrt(var_g),
    log_marginal_likelihood = lml
  )
  colnames(values) <- c(
    "lp__", sprintf("alpha[%d]", seq_len(data$M)), "lambda",
    sprintf("log_g[%d]", seq_len(data$L)), "mean_g", "var_g", "sd_g",
    "log_marginal_likelihood"
  )
  # Iterations by chains by variables: the rows of `values` are the draws of
  # the first chain and then those of the second.
  array_draws <- posterior::as_draws_array(
    array(values, dim = c(n_iter, n_chains, ncol(values)),
          dimnames = list(NULL, NULL, colnames(values)))
  )
  list(
    draws = function(variables = NULL, format = "draws_array") array_draws,
    values = values
  )
}

programs_case <- function(tau_c = 0, n_draws = 40L) {
  withr::local_seed(4046)
  theta_hat <- c(-0.45, -0.1, 0, 0.25, 0.55, 0.9)
  sigma <- c(0.12, 0.18, 0.15, 0.22, 0.2, 0.3)
  grid <- make_efron_grid(theta_hat, sigma, L = 51L, M = 4L, grid_method = "paper_realdata")
  data <- list(K = length(theta_hat), theta_hat = theta_hat, sigma = sigma, L = grid$L,
               grid = grid$grid, M = grid$M, B = unname(grid$B), tau_c = tau_c,
               store_grid_quantities = 0L)
  alpha <- matrix(stats::rnorm(n_draws * grid$M, sd = 0.8), nrow = n_draws)
  list(data = data, fit = programs_fake_fit(alpha, data))
}

test_that("the completed draws have the variables of version 0.3.0 in its order", {
  case <- programs_case()
  K <- case$data$K
  L <- case$data$L
  M <- case$data$M
  compact <- programs_ns(".bef_complete_mix_draws")(case$fit$draws(), case$data)
  expect_s3_class(compact, "draws_array")
  expect_identical(dim(compact)[1:2], c(20L, 2L))
  expect_identical(
    dimnames(compact)[[3L]],
    c("lp__", sprintf("alpha[%d]", seq_len(M)), "lambda", sprintf("log_g[%d]", seq_len(L)),
      "mean_g", "var_g", "sd_g", "effective_params", "log_marginal_likelihood")
  )
  full <- programs_ns(".bef_complete_mix_draws")(case$fit$draws(), case$data,
                                                site_draws = TRUE, seed = 5L)
  expect_identical(
    dimnames(full)[[3L]],
    c("lp__", sprintf("alpha[%d]", seq_len(M)), "lambda", sprintf("log_g[%d]", seq_len(L)),
      "mean_g", "var_g", "sd_g",
      as.vector(vapply(c("theta_map", "theta_mean", "theta_sd", "theta_rep", "log_lik"),
                       function(name) sprintf("%s[%d]", name, seq_len(K)), character(K))),
      "effective_params", "log_marginal_likelihood")
  )
  # What CmdStan wrote is carried over unchanged, draw by draw.
  written <- colnames(case$fit$values)
  expect_identical(programs_plain(full)[, written], programs_plain(case$fit$draws())[, written])
  expect_identical(programs_plain(full)[, written], unname(case$fit$values)[, seq_along(written)] |>
                     structure(dimnames = list(NULL, written)))
  expect_identical(programs_plain(compact), programs_plain(full)[, dimnames(compact)[[3L]]])
})

test_that("the site quantities added in R belong to the draw they are stored with", {
  for (tau_c in c(0, 0.05)) {
    case <- programs_case(tau_c = tau_c)
    data <- case$data
    full <- programs_plain(
      programs_ns(".bef_complete_mix_draws")(case$fit$draws(), data, site_draws = TRUE,
                                            seed = 5L)
    )
    shrink <- tau_c^2 / (tau_c^2 + data$sigma^2)
    for (row in c(1L, 7L, 20L, 21L, 40L)) {
      log_g <- full[row, sprintf("log_g[%d]", seq_len(data$L))]
      mean_i <- sd_i <- lik_i <- map_i <- numeric(data$K)
      for (i in seq_len(data$K)) {
        log_post <- log_g + stats::dnorm(data$theta_hat[[i]], data$grid,
                                         sqrt(data$sigma[[i]]^2 + tau_c^2), log = TRUE)
        w <- exp(log_post - max(log_post))
        lik_i[[i]] <- max(log_post) + log(sum(w))
        w <- w / sum(w)
        centers <- data$grid + shrink[[i]] * (data$theta_hat[[i]] - data$grid)
        mean_i[[i]] <- sum(w * centers)
        sd_i[[i]] <- sqrt(shrink[[i]] * data$sigma[[i]]^2 + sum(w * (centers - mean_i[[i]])^2))
        map_i[[i]] <- centers[[which.max(log_post)]]
      }
      at <- function(name) unname(full[row, sprintf("%s[%d]", name, seq_len(data$K))])
      expect_equal(at("theta_mean"), mean_i, tolerance = 1e-10)
      expect_equal(at("theta_sd"), sd_i, tolerance = 1e-9)
      expect_equal(at("theta_map"), map_i, tolerance = 1e-12)
      expect_equal(at("log_lik"), lik_i, tolerance = 1e-10)
      expect_equal(unname(full[row, "effective_params"]), sum(sd_i^2 / data$sigma^2),
                   tolerance = 1e-9)
      expect_equal(unname(full[row, "log_marginal_likelihood"]), sum(lik_i), tolerance = 1e-10)
    }
  }
})

test_that("the draws of the effects are reproducible and leave R's generator alone", {
  case <- programs_case()
  complete <- function(seed) {
    programs_plain(programs_ns(".bef_complete_mix_draws")(
      case$fit$draws(), case$data, site_draws = TRUE, seed = seed
    ))
  }
  withr::local_seed(99)
  before <- .Random.seed
  first <- complete(5L)
  expect_identical(.Random.seed, before)
  expect_identical(complete(5L), first)
  other <- complete(6L)
  rep_cols <- sprintf("theta_rep[%d]", seq_len(case$data$K))
  expect_false(identical(other[, rep_cols], first[, rep_cols]))
  expect_identical(other[, setdiff(colnames(first), rep_cols)],
                   first[, setdiff(colnames(first), rep_cols)])
  # With point masses every draw of an effect is a grid point.
  expect_true(all(first[, rep_cols] %in% case$data$grid))

  # The draws follow the posterior probabilities of the grid points.
  big <- programs_case(n_draws = 4000L)
  alpha <- matrix(rep(c(0.3, -0.4, 0.2, 0.1), each = 4000L), nrow = 4000L)
  fit <- programs_fake_fit(alpha, big$data)
  draws <- programs_plain(programs_ns(".bef_complete_mix_draws")(
    fit$draws(), big$data, site_draws = TRUE, seed = 11L
  ))
  for (i in c(1L, 4L)) {
    rep_i <- draws[, sprintf("theta_rep[%d]", i)]
    expect_lt(abs(mean(rep_i) - draws[1L, sprintf("theta_mean[%d]", i)]),
              4 * draws[1L, sprintf("theta_sd[%d]", i)] / sqrt(4000))
    expect_equal(stats::sd(rep_i), unname(draws[1L, sprintf("theta_sd[%d]", i)]),
                 tolerance = 0.06)
  }
})

test_that("post-processing gives the same fit pieces for the two kinds of storage", {
  case <- programs_case()
  post <- function(site_draws) {
    suppressWarnings(postprocess_stan_draws(
      case$fit, case$data, site_draws = site_draws, program = "efron_mix", seed = 3L
    ))
  }
  compact <- post(FALSE)
  full <- post(TRUE)
  expect_identical(compact$posterior, full$posterior)
  expect_identical(compact$theta_summary, full$theta_summary)
  expect_named(compact$posterior, c("mean_g", "var_g", "sd_g", "effective_params",
                                    "log_marginal_likelihood"))
  expect_false(any(grepl("^theta_", dimnames(compact$draws)[[3L]])))
  expect_true(all(is.finite(compact$posterior$effective_params)))
  expect_true(all(compact$posterior$effective_params > 0))
})

test_that("with point masses the two programs give the same posterior for the callback data", {
  skip_if_no_cmdstan()
  run <- function(settings) {
    context <- programs_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)),
      theta_hat = krw_discrimination$yi, sigma = krw_discrimination$sei,
      grid_method = "paper_realdata", iter_sampling = 1000L, seed = 20251L,
      settings = settings
    )
    suppressWarnings(programs_ns(".bef_fit_assemble")(programs_ns(".bef_fit_run_sampler")(context)))
  }
  old <- run("0.3.0")
  new <- run(NULL)
  expect_identical(attr(old$metadata, "program", exact = TRUE), "efron_re")
  expect_identical(attr(new$metadata, "program", exact = TRUE), "efron_mix")
  expect_false(identical(old$metadata$stan_file_sha256, new$metadata$stan_file_sha256))
  expect_identical(dimnames(new$draws)[[3L]], dimnames(old$draws)[[3L]])
  expect_identical(attr(new$metadata, "sampler_settings", exact = TRUE)$sig_figs, 12L)
  expect_null(attr(old$metadata, "sampler_settings", exact = TRUE)$sig_figs)

  # Each quantity agrees within Monte Carlo error. Four standard errors and
  # not two: the chains differ between platforms and CmdStan versions, and
  # a test at two would fail by chance once in twenty runs.
  for (name in c("mean_g", "sd_g", "effective_params", "log_marginal_likelihood")) {
    a <- posterior::subset_draws(old$draws, variable = name)
    b <- posterior::subset_draws(new$draws, variable = name)
    se <- sqrt(posterior::mcse_mean(a)^2 + posterior::mcse_mean(b)^2)
    expect_lt(abs(mean(a) - mean(b)), 4 * se)
  }
  # So does the posterior mean of every site: five standard errors, for 97
  # comparisons.
  site_means <- function(fit) {
    conditional <- programs_ns(".bef_site_conditional")(programs_ns(".bef_engine")(fit))$mean
    chains <- dim(fit$draws)[[2L]]
    list(
      mean = colMeans(conditional),
      mcse = apply(conditional, 2L, function(x) {
        posterior::mcse_mean(matrix(x, ncol = chains))
      })
    )
  }
  a <- site_means(old)
  b <- site_means(new)
  expect_equal(a$mean, old$metadata$theta_summary$mean, tolerance = 1e-10)
  expect_true(all(abs(a$mean - b$mean) < 5 * sqrt(a$mcse^2 + b$mcse^2)))

  # Compact and full storage of the new program share the sampled draws.
  full <- {
    context <- programs_ns(".bef_fit_prepare_context")(
      call = quote(bayes_efron_fit(theta_hat, sigma)),
      theta_hat = krw_discrimination$yi, sigma = krw_discrimination$sei,
      grid_method = "paper_realdata", iter_sampling = 1000L, seed = 20251L,
      site_draws = TRUE
    )
    suppressWarnings(programs_ns(".bef_fit_assemble")(programs_ns(".bef_fit_run_sampler")(context)))
  }
  shared <- dimnames(new$draws)[[3L]]
  expect_identical(programs_plain(full$draws)[, shared], programs_plain(new$draws))
  # And so every summary of the sites and of g.
  expect_identical(site_summary(full), site_summary(new))
  expect_identical(g_summary(full), g_summary(new))
  expect_identical(coef(full), coef(new))
  expect_equal(confint(full), confint(new))
  expect_equal(dim(full$draws)[[3L]], length(shared) + 5L * nrow(krw_discrimination))
})
