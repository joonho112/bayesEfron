// ==============================================================================
// Fully Bayesian Meta-Analytic Deconvolution with Efron's Log-Spline Prior
//
// Reference: Lee & Sui (2025). "Fully Bayesian Inference for Meta-Analytic 
//            Deconvolution Using Efron's Log-Spline Prior"
//
// Model Description:
// This Stan program implements a fully Bayesian approach to meta-analytic
// deconvolution using Efron's log-spline prior. The model estimates a smooth
// prior distribution g(θ) from noisy observations while properly accounting
// for both sampling uncertainty and prior uncertainty.
//
// Key Features:
// - Flexible nonparametric prior via log-spline representation
// - Automatic regularization through hierarchical Bayesian framework
// - Proper uncertainty quantification for site-specific effects
// - Numerical stability through log-space computations
//
// File: efron_re.stan
// Package release: bayesEfron 0.1.0
// Model family: RE
// ==============================================================================

data {
  // Site-level data
  int<lower=1> K;                  // Number of sites/studies
  vector[K] theta_hat;              // Observed effect estimates
  vector<lower=0>[K] sigma;         // Standard errors of estimates
  
  // Discretization grid for prior g
  int<lower=1> L;                  // Number of grid points (typically 101)
  vector[L] grid;                   // Grid points spanning support of θ
  
  // Spline basis specification
  int<lower=1> M;                  // Degrees of freedom for splines (typically 6)
  matrix[L, M] B;                   // Natural cubic spline basis matrix
}

parameters {
  // Spline coefficients for log-density representation
  vector[M] alpha;                  // Coefficients for log g(θ)
  
  // Regularization parameter
  real<lower=0> lambda;             // Controls smoothness of prior g
}

transformed parameters {
  // Log-spline representation of prior
  vector[L] log_w = B * alpha;      // Log unnormalized density at grid points
  
  // Normalized prior distribution
  vector[L] log_g = log_softmax(log_w);  // Log normalized density (sums to 1)
  simplex[L] g = softmax(log_w);         // Normalized density on simplex
}

model {
  // ============================================================================
  // PRIORS
  // ============================================================================
  
  // Hyperprior on regularization parameter
  // Half-Cauchy(0, 5) is weakly informative, allowing adaptation to data
  lambda ~ cauchy(0, 5);
  
  // Conditional prior on spline coefficients
  // Ridge penalty with variance 1/lambda encourages smoothness
  alpha ~ normal(0, inv_sqrt(lambda));
  
  // ============================================================================
  // LIKELIHOOD
  // ============================================================================
  
  // Mixture likelihood for observed effects
  // Each θ_hat_i comes from mixture: Σ_j g_j × Normal(grid_j, sigma_i)
  for (i in 1:K) {
    vector[L] log_components;
    
    // Compute log-likelihood for each mixture component
    for (j in 1:L) {
      log_components[j] = log_g[j]  // Prior weight for component j
                        + normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
    }
    
    // Add log marginal likelihood using log-sum-exp for numerical stability
    target += log_sum_exp(log_components);
  }
}

generated quantities {
  // ============================================================================
  // PRIOR DISTRIBUTION SUMMARIES
  // ============================================================================
  
  // Moments of estimated prior g(θ)
  real mean_g = dot_product(g, grid);                     // E[θ]
  real var_g = dot_product(g, square(grid - mean_g));     // Var[θ]
  real sd_g = sqrt(var_g);                                // SD[θ]
  
  // ============================================================================
  // POSTERIOR SUMMARIES FOR SITE-SPECIFIC EFFECTS
  // ============================================================================
  
  // Three types of posterior summaries for each site
  vector[K] theta_map;   // Maximum a posteriori (MAP) estimates
  vector[K] theta_mean;  // Posterior means (optimal under squared error loss)
  vector[K] theta_rep;   // Posterior draws (for credible intervals)
  
  // Additional posterior uncertainty measures
  vector[K] theta_sd;    // Posterior standard deviations
  
  for (i in 1:K) {
    // ------------------------------------------------------------------------
    // Compute posterior distribution for site i
    // ------------------------------------------------------------------------
    
    vector[L] log_post;  // Log posterior at each grid point
    
    // Bayes' rule: posterior ∝ prior × likelihood
    for (j in 1:L) {
      log_post[j] = log_g[j]  // Log prior
                  + normal_lpdf(theta_hat[i] | grid[j], sigma[i]);  // Log likelihood
    }
    
    // ------------------------------------------------------------------------
    // MAP estimate (posterior mode)
    // ------------------------------------------------------------------------
    
    int max_idx = 1;
    for (j in 2:L) {
      if (log_post[j] > log_post[max_idx]) {
        max_idx = j;
      }
    }
    theta_map[i] = grid[max_idx];
    
    // ------------------------------------------------------------------------
    // Normalize posterior to get weights
    // ------------------------------------------------------------------------
    
    // Use log-sum-exp trick for numerical stability
    real log_post_max = max(log_post);
    vector[L] w = exp(log_post - log_post_max);
    w = w / sum(w);  // Now w contains posterior probabilities
    
    // ------------------------------------------------------------------------
    // Posterior mean
    // ------------------------------------------------------------------------
    
    theta_mean[i] = dot_product(w, grid);
    
    // ------------------------------------------------------------------------
    // Posterior standard deviation
    // ------------------------------------------------------------------------
    
    real second_moment = dot_product(w, square(grid));
    theta_sd[i] = sqrt(second_moment - square(theta_mean[i]));
    
    // ------------------------------------------------------------------------
    // Posterior draw (for uncertainty quantification)
    // ------------------------------------------------------------------------
    
    // Sample from discrete posterior distribution
    theta_rep[i] = grid[categorical_rng(w)];
  }
  
  // ============================================================================
  // MODEL DIAGNOSTICS AND CHECKS
  // ============================================================================
  
  // Effective number of parameters (for model comparison)
  real effective_params;
  {
    vector[K] posterior_vars;
    for (i in 1:K) {
      // Compute posterior variance for each site
      vector[L] log_post;
      for (j in 1:L) {
        log_post[j] = log_g[j] + normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
      }
      real log_post_max = max(log_post);
      vector[L] w = exp(log_post - log_post_max);
      w = w / sum(w);
      
      real post_mean = dot_product(w, grid);
      real post_second_moment = dot_product(w, square(grid));
      posterior_vars[i] = post_second_moment - square(post_mean);
    }
    
    // Effective parameters = K - sum of shrinkage factors
    effective_params = K - sum(posterior_vars ./ square(sigma));
  }
  
  // Log marginal likelihood (for model comparison)
  real log_marginal_likelihood = 0;
  for (i in 1:K) {
    vector[L] log_components;
    for (j in 1:L) {
      log_components[j] = log_g[j] + normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
    }
    log_marginal_likelihood += log_sum_exp(log_components);
  }
}
