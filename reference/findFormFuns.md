# `findFormFuns` used by [averageObs](averageObs.md) to calculate proper averages

The purpose is to properly derive data for the average observation in
the data by being 'aware' of formulas that contain interactions and/or
function calls. For example, in the old behavior, if the formula
contained a square term specified as `I(x^2)`, we were returning the
mean of `x(^2)` not the square of mean(x).

## Usage

``` r
findFormFuns(merMod, origData = NULL)
```

## Arguments

- merMod:

  the merMod object from which to draw the average observation

- origData:

  (default=NULL) a data frame containing the original, untransformed
  data used to call the model. This MUST be specified if the original
  variables used in formula function calls are NOT present as 'main
  effects'.

## Value

a data frame with a single row for the average observation, but with
full factor levels. See details for more.

## Details

Matrix-valued response columns (e.g. the `cbind(successes, failures)`
left-hand side of a binomial GLMM) are detected and dropped from the
working frame before averaging, since they cannot be collapsed to a
single scalar. The returned frame therefore has no response column for
matrix-LHS models; see [`averageObs`](averageObs.md) for rationale.
