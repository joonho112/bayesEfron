// Fully Bayesian deconvolution with Efron's log-spline prior
//
// Lee and Sui (2025), "Fully Bayesian Inference for Meta-Analytic
// Deconvolution Using Efron's Log-Spline Prior", Mathematics 13(16), 2639.
//
// For sites i = 1, ..., K (the studies of a meta-analysis, the sites of a
// multisite trial),
//
//   theta_hat[i] | theta[i] ~ normal(theta[i], sigma[i]),   theta[i] ~ g,
//
// where g is a distribution on a grid of L points with
//
//   log g = log_softmax(B * alpha),
//
// B is a natural cubic spline basis with M columns, and
//
//   alpha[m] | lambda ~ normal(0, 1 / sqrt(lambda)),   lambda ~ half-Cauchy(0, 5).
//
// The theta[i] are summed out on the grid, so alpha and lambda are the only
// sampled parameters.

data {
  // Site-level data
  int<lower=1> K;                  // number of sites
  vector[K] theta_hat;              // effect estimates
  vector<lower=0>[K] sigma;         // standard errors of the estimates
  
  // Grid on which g is defined
  int<lower=1> L;                  // number of grid points
  vector[L] grid;                   // grid points
  
  // Spline basis
  int<lower=1> M;                  // number of basis functions
  matrix[L, M] B;                   // natural cubic spline basis on the grid

  // 1 also writes log_w = B * alpha and g = exp(log_g) to the output. Both
  // follow from alpha and log_g, so the default, 0, leaves them out.
  int<lower=0, upper=1> store_grid_quantities;
}

transformed data {
  // normal_lpdf(theta_hat[i] | grid[j], sigma[i]) depends on the data only, so
  // it is computed once here and not at every gradient evaluation.
  matrix[K, L] log_lik_grid;
  for (i in 1:K) {
    for (j in 1:L) {
      log_lik_grid[i, j] = normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
    }
  }

  // Shifted likelihood for the marginal in the model block:
  //
  //   log_sum_exp_j( LL[i,j] + log_g[j] )
  //     = row_max[i] + log( sum_j exp(LL[i,j] - row_max[i]) * g[j] )
  //     = row_max[i] + log( lik_shifted[i,] * g ).
  //
  // row_max depends on the data only, so the marginal likelihood of all the
  // sites is one matrix-vector product. Every row of lik_shifted has maximum
  // 1, so the product for a site is at least g at the grid point nearest its
  // estimate. It rounds to zero only if g underflows there, which needs
  // coefficients far outside the range that the sampler visits. The log
  // density is then minus infinity, and Stan does not accept such a point.
  vector[K] log_lik_row_max;
  matrix[K, L] lik_grid_shifted;
  for (i in 1:K) {
    log_lik_row_max[i] = max(log_lik_grid[i]);
    lik_grid_shifted[i] = exp(log_lik_grid[i] - log_lik_row_max[i]);
  }
  real log_lik_row_max_sum = sum(log_lik_row_max);
}

parameters {
  // Spline coefficients of log g
  vector[M] alpha;                  // coefficients of the basis functions
  
  // Prior precision of alpha; larger values pull g toward the uniform
  // distribution on the grid
  real<lower=0> lambda;             // precision
}

transformed parameters {
  // g on the log scale. Only log_g is stored; g = exp(log_g) is formed where
  // it is needed.
  vector[L] log_g = log_softmax(B * alpha);
}

model {
  // Priors
  
  // Half-Cauchy(0, 5), because lambda is constrained to be positive
  lambda ~ cauchy(0, 5);
  
  // Normal with variance 1 / lambda
  alpha ~ normal(0, inv_sqrt(lambda));
  
  // Likelihood
  
  // Marginal likelihood of theta_hat with theta summed out on the grid:
  // theta_hat[i] has density sum_j g[j] * normal(theta_hat[i] | grid[j], sigma[i]).
  // See the transformed data block for the identity used here.
  target += log_lik_row_max_sum + sum(log(lik_grid_shifted * exp(log_g)));
}

generated quantities {
  // Everything below is computed in one pass over the sites. Only the
  // quantities declared at this level are written to the output; working
  // values are local to the block at the bottom.

  // Mean, variance and standard deviation of g
  real mean_g;
  real var_g;
  real sd_g;

  // Posterior of theta[i] given alpha
  vector[K] theta_map;   // posterior mode on the grid
  vector[K] theta_mean;  // posterior mean
  vector[K] theta_sd;    // posterior standard deviation
  vector[K] theta_rep;   // one draw from the posterior
  vector[K] log_lik;     // marginal log density of theta_hat[i]

  // Effective number of parameters (defined below)
  real effective_params;

  // Sum of log_lik
  real log_marginal_likelihood;

  // Written only when store_grid_quantities is 1
  vector[store_grid_quantities ? L : 0] log_w;
  vector[store_grid_quantities ? L : 0] g;

  {
    vector[L] g_local = exp(log_g);
    vector[K] posterior_vars;

    mean_g = dot_product(g_local, grid);
    var_g = dot_product(g_local, square(grid - mean_g));
    sd_g = sqrt(var_g);

    for (i in 1:K) {
      // Bayes' rule on the grid, computed once and reused for every summary
      vector[L] log_post = to_vector(log_lik_grid[i]) + log_g;
      real log_post_max = max(log_post);
      vector[L] w_unnorm = exp(log_post - log_post_max);
      real w_sum = sum(w_unnorm);
      vector[L] w = w_unnorm / w_sum;

      // Marginal log density: log_sum_exp(log_post), from the pieces above
      log_lik[i] = log_post_max + log(w_sum);

      // Posterior mode
      int max_idx = 1;
      for (j in 2:L) {
        if (log_post[j] > log_post[max_idx]) {
          max_idx = j;
        }
      }
      theta_map[i] = grid[max_idx];

      // Posterior mean and variance. The variance is taken about the mean,
      // and not as a difference of two moments, which loses its digits when
      // the effects are far from zero.
      real m1 = dot_product(w, grid);
      theta_mean[i] = m1;
      posterior_vars[i] = dot_product(w, square(grid - m1));
      theta_sd[i] = sqrt(posterior_vars[i]);

      // Posterior draw on the grid
      theta_rep[i] = grid[categorical_rng(w)];
    }

    // Effective number of parameters. In the normal-normal model the posterior
    // mean of theta[i] is w[i] * theta_hat[i] + (1 - w[i]) * mu with
    // w[i] = tau^2 / (sigma[i]^2 + tau^2), and the posterior variance is
    // V[i] = sigma[i]^2 * w[i]. So V[i] / sigma[i]^2 = w[i] is the weight the
    // data carry for site i, and the sum over sites is the trace of the
    // smoother, which grows as the estimates are shrunk less. The same ratio
    // is used here with the posterior variance on the grid. In the
    // normal-normal model it lies between 0 and K; with a general g a single
    // term can exceed 1.
    effective_params = sum(posterior_vars ./ square(sigma));

    log_marginal_likelihood = sum(log_lik);

    if (store_grid_quantities) {
      log_w = B * alpha;
      g = g_local;
    }
  }
}
