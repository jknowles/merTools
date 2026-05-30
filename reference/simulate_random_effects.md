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
  seed = NULL,
  new.levels = c("zero", "draw")
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

- new.levels:

  Character; how to treat grouping levels present in \`newdata\` but
  absent from the fitted model. \`"zero"\` drops the random effect for
  such levels; \`"draw"\` samples each unobserved level's effect from
  the estimated random-effect covariance (\`VarCorr\`).

## Value

List of matrices (rows = observations, cols = simulations).
