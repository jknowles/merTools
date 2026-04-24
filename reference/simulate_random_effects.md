# Simulate random‑effect contributions for all grouping factors

For each random effect term the function draws \`n.sims\` samples from
the conditional multivariate normal distribution (BLUPs + post‑var). The
result is a named list where each element is an array of dimensions
\`nrow(newdata)\` × \`n.sims\`.

## Usage

``` r
simulate_random_effects(
  merMod,
  newdata,
  n.sims,
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

- .parallel:

  Logical flag to enable parallel execution via foreach.

- seed:

  Optional random seed for reproducibility.

## Value

List of matrices (rows = observations, cols = simulations).
