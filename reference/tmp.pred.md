# Compute random effect predictions for a single grouping factor

Internal helper that multiplies predictor data by simulated random
effect coefficients for each grouping level. Handles new levels not
present in the original model data by returning zeros for their random
effect contributions.

## Usage

``` r
tmp.pred(data, coefs, group)
```

## Arguments

- data:

  data.frame with predictor columns and a grouping factor column. The
  last column must be the grouping factor.

- coefs:

  3D array of simulated coefficients with dimensions (n.sims x terms x
  levels).

- group:

  character name of the grouping factor column in data.

## Value

matrix of predicted values with dimensions (nrow(data) x n.sims).
