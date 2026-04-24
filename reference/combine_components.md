# Combine fixed, random and residual components into a prediction array

Combine fixed, random and residual components into a prediction array

## Usage

``` r
combine_components(
  fixed_mat,
  random_list,
  sigma_vec,
  include.resid.var = TRUE,
  which = c("full", "fixed", "random", "all"),
  family = NULL,
  link = NULL,
  weights = NULL,
  use.probability = FALSE
)
```

## Arguments

- fixed_mat:

  Matrix of dimensions \`nrow(newdata)\` × \`n.sims\`.

- random_list:

  List returned by \`simulate_random_effects()\`.

- sigma_vec:

  Residual standard deviations (output of
  \`simulate_residual_variance\`).

- include.resid.var:

  Logical, whether to add residual variance.

- which:

  Character: "full", "fixed", "random" or "all".

## Value

Array \`nrow(newdata)\` × \`n.sims\` (or list when \`which="all"\`).
