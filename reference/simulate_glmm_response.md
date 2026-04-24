# Simulate GLMM responses from conditional distribution

Internal helper that simulates responses from the conditional
distribution of Y\|b (given random effects) for GLMMs. This implements
the correct theoretical approach where residual variance is inherent in
the distribution.

## Usage

``` r
simulate_glmm_response(
  yhat,
  family,
  link,
  weights = NULL,
  use.probability = FALSE
)
```

## Arguments

- yhat:

  Matrix of linear predictors (fixed + random effects).

- family:

  Character string specifying the distribution family.

- link:

  Character string specifying the link function.

- weights:

  Optional weights (for binomial with n-trial parameter).

- use.probability:

  For binomial, return proportions (0-1) instead of counts?

## Value

Matrix of simulated responses with same dimensions as yhat.
