// Fully Bayesian deconvolution with Efron's log-spline prior
//
// Lee and Sui (2025), "Fully Bayesian Inference for Meta-Analytic
// Deconvolution Using Efron's Log-Spline Prior", Mathematics 13(16), 2639.
//
// This is the program of version 0.4. The model is that of efron_re.stan,
// the program of versions up to 0.3, with one generalization: g has a
// component at each grid point that is either a point mass (tau_c = 0, the
// model of the paper) or a normal distribution with standard deviation
// tau_c. For sites i = 1, ..., K,
//
//   theta_hat[i] | theta[i] ~ normal(theta[i], sigma[i]),   theta[i] ~ g,
//
//   g = sum over j of g[j] * normal(grid[j], tau_c),   log g = log_softmax(B * alpha),
//
// B is a natural cubic spline basis with M columns, and
//
//   alpha[m] | lambda ~ normal(0, 1 / sqrt(lambda)),   lambda ~ half-Cauchy(0, 5).
//
// The theta[i] are integrated out: given component j, theta_hat[i] is normal
// with mean grid[j] and variance sigma[i]^2 + tau_c^2. So alpha and lambda
// are the only sampled parameters.
//
// Unlike efron_re.stan, this program writes nothing for the sites. Every
// summary of a site follows from log_g and the data, and the package
// computes it in R from the draws of log_g.

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

  // Standard deviation of the component at each grid point; 0 for point
  // masses
  real<lower=0> tau_c;

  // 1 also writes log_w = B * alpha and g = exp(log_g) to the output. Both
  // follow from alpha and log_g, so the default, 0, leaves them out.
  int<lower=0, upper=1> store_grid_quantities;
}

transformed data {
  // Shifted likelihood for the marginal in the model block. With
  // LL[i,j] = normal_lpdf(theta_hat[i] | grid[j], sqrt(sigma[i]^2 + tau_c^2)),
  //
  //   log_sum_exp_j( LL[i,j] + log_g[j] )
  //     = row_max[i] + log( sum_j exp(LL[i,j] - row_max[i]) * g[j] )
  //     = row_max[i] + log( lik_shifted[i,] * g ).
  //
  // LL depends on the data only, so it is computed once here, and the
  // marginal likelihood of all the sites is one matrix-vector product. Every
  // row of lik_shifted has maximum 1, so the product for a site is at least g
  // at the grid point nearest its estimate. It rounds to zero only if g
  // underflows there, which needs coefficients far outside the range that
  // the sampler visits. The log density is then minus infinity, and Stan
  // does not accept such a point.
  matrix[K, L] lik_grid_shifted;
  real log_lik_row_max_sum;
  {
    vector[K] sd_total = sqrt(square(sigma) + square(tau_c));
    vector[K] log_lik_row_max;
    for (i in 1:K) {
      row_vector[L] log_lik_row;
      for (j in 1:L) {
        log_lik_row[j] = normal_lpdf(theta_hat[i] | grid[j], sd_total[i]);
      }
      log_lik_row_max[i] = max(log_lik_row);
      lik_grid_shifted[i] = exp(log_lik_row - log_lik_row_max[i]);
    }
    log_lik_row_max_sum = sum(log_lik_row_max);
  }
}

parameters {
  // Spline coefficients of log g
  vector[M] alpha;                  // coefficients of the basis functions

  // Prior precision of alpha; larger values pull g toward the uniform
  // distribution on the grid
  real<lower=0> lambda;             // precision
}

transformed parameters {
  // The probabilities of the components on the log scale. Only log_g is
  // stored; g = exp(log_g) is formed where it is needed.
  vector[L] log_g = log_softmax(B * alpha);
}

model {
  // Priors

  // Half-Cauchy(0, 5), because lambda is constrained to be positive
  lambda ~ cauchy(0, 5);

  // Normal with variance 1 / lambda
  alpha ~ normal(0, inv_sqrt(lambda));

  // Likelihood

  // Marginal likelihood of theta_hat with theta integrated out. See the
  // transformed data block for the identity used here.
  target += log_lik_row_max_sum + sum(log(lik_grid_shifted * exp(log_g)));
}

generated quantities {
  // Mean, variance and standard deviation of g. A component adds its own
  // variance to that of the grid points.
  real mean_g;
  real var_g;
  real sd_g;

  // Marginal log density of the estimates given alpha
  real log_marginal_likelihood;

  // Written only when store_grid_quantities is 1
  vector[store_grid_quantities ? L : 0] log_w;
  vector[store_grid_quantities ? L : 0] g;

  {
    vector[L] g_local = exp(log_g);

    mean_g = dot_product(g_local, grid);
    var_g = dot_product(g_local, square(grid - mean_g)) + square(tau_c);
    sd_g = sqrt(var_g);

    log_marginal_likelihood =
      log_lik_row_max_sum + sum(log(lik_grid_shifted * g_local));

    if (store_grid_quantities) {
      log_w = B * alpha;
      g = g_local;
    }
  }
}
