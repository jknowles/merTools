# Simulate residual variance for merMod objects

Internal helper that reproduces the sigma-drawing logic from the
original \`predictInterval()\` implementation. It returns a numeric
vector of length \`n.sims\` containing the simulated residual standard
deviations.

## Usage

``` r
simulate_residual_variance(merMod, n.sims)
```

## Arguments

- merMod:

  A merMod object from lme4.

- n.sims:

  Number of simulation draws.

## Value

Numeric vector of length \`n.sims\`.
