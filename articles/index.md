# Articles

### Case studies

- [Twin Towers: recovering a distribution hidden by
  noise](https://joonho112.github.io/bayesEfron/articles/twin-towers.md):

  A simulated data set of 1,500 sites with known true effects in two
  separated groups, fitted from the estimates and standard errors alone.

- [Nineteen studies of teacher
  expectancy](https://joonho112.github.io/bayesEfron/articles/teacher-expectancy.md):

  A case study of a meta-analysis with 19 studies: the estimated
  distribution of effects, its sensitivity to the range of the grid, and
  a simulation study of fits of this size.

- [Small classes in 79
  schools](https://joonho112.github.io/bayesEfron/articles/class-size.md):

  Analyze a multisite trial in education: build an estimate and a
  standard error for each of 79 schools from student data, fit the
  model, and use alternative standard errors and a simulation with known
  effects to judge the evidence that the effect of small classes varies
  across schools.

- [bayesEfron and empirical Bayes methods
  compared](https://joonho112.github.io/bayesEfron/articles/empirical-bayes.md):

  bayesEfron and the empirical Bayes methods in the ebnm package, fitted
  to the same data: agreement of the posterior means, differences in the
  intervals and in the estimated distribution of effects, and the
  contribution of the uncertainty about that distribution.

### Guides

- [Estimates, standard errors and site
  labels](https://joonho112.github.io/bayesEfron/articles/data-input.md):

  Getting estimates and standard errors into a fit from a data frame,
  two vectors or metafor, and keeping the names of the sites through to
  the tables and plots.

- [Plots of a fitted
  model](https://joonho112.github.io/bayesEfron/articles/plots.md):

  A gallery of the plot functions: the question that each answers, the
  call, and what the result shows.

- [Choosing the grid and checking sensitivity to
  it](https://joonho112.github.io/bayesEfron/articles/choosing-a-grid.md):

  The choices that define the grid, the default method for its range,
  checks of the range after a fit, and a comparison of results across
  grids.

- [Checking a
  fit](https://joonho112.github.io/bayesEfron/articles/diagnostics.md):

  Sampler diagnostics, grid diagnostics, a posterior predictive check,
  and how to respond to a warning.

- [The model and the posterior distribution of a site
  effect](https://joonho112.github.io/bayesEfron/articles/model.md):

  The deconvolution problem, the log-spline family for the distribution
  of effects, its prior, and the posterior distribution of a site
  effect.

- [Reproducibility and interval
  coverage](https://joonho112.github.io/bayesEfron/articles/implementation.md):

  Reproducing a fit, the draws stored in a fit, the two Stan programs,
  checks of the computations, interval coverage in simulated data, and
  the relation to the paper.
