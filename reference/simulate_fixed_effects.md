# Simulate fixed-effect predictions

Draws \`n.sims\` samples from the multivariate normal distribution
defined by the fixed-effect estimates and their variance-covariance
matrix, then multiplies by the model matrix to produce predictions. It
also respects \`ignore.fixed.terms\` and optional intercept variance
adjustment.

## Usage

``` r
simulate_fixed_effects(
  merMod,
  newdata,
  n.sims,
  ignore.fixed.terms = NULL,
  fix.intercept.variance = FALSE,
  .parallel = FALSE,
  seed = NULL
)
```

## Arguments

- merMod:

  A merMod object.

- newdata:

  Data frame containing the prediction covariates.

- n.sims:

  Number of simulations.

- ignore.fixed.terms:

  Vector of fixed-effect names or indices to treat as known.

- fix.intercept.variance:

  Logical flag for intercept variance correction.

- .parallel:

  Logical flag to enable parallel execution via foreach.

## Value

Matrix of dimensions \`nrow(newdata)\` x \`n.sims\` containing simulated
fixed-effect predictions.
