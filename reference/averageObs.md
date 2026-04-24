# Find the average observation for a merMod object

Extract a data frame of a single row that represents the average
observation in a merMod object. This function also allows the user to
pass a series of conditioning argument to calculate the average
observation conditional on other characteristics.

## Usage

``` r
averageObs(merMod, varList = NULL, origData = NULL, ...)
```

## Arguments

- merMod:

  a merMod object

- varList:

  optional, a named list of conditions to subset the data on

- origData:

  (default=NULL) a data frame containing the original, untransformed
  data used to call the model. This MUST be specified if the original
  variables used in formula function calls are NOT present as 'main
  effects'.

- ...:

  not used currently

## Value

a data frame with a single row for the average observation, but with
full factor levels. See details for more.

## Details

Each character and factor variable in the data.frame is assigned to the
modal category and each numeric variable is collapsed to the mean.
Currently if mode is a tie, returns a "." Uses the collapseFrame
function.

For models with a scalar left-hand side (e.g. `lmer(y ~ ...)`), the
response column is included in the output and is set to the mean of the
observed response. For models with a matrix-valued left-hand side – most
commonly two-column [`cbind()`](https://rdrr.io/r/base/cbind.html)
specifications in binomial GLMMs such as
`glmer(cbind(successes, failures) ~ ..., family = binomial)` – the
response column is omitted from the output. A matrix response cannot be
meaningfully collapsed to a single "average" value, and `averageObs()`
is primarily intended to produce `newdata` for
[`predict`](https://rdrr.io/r/stats/predict.html) /
[`predictInterval`](predictInterval.md), both of which ignore the
response column in `newdata`. Callers that iterate over
`names(averageObs(merMod))` or compare against `merMod@frame` should not
assume column parity for matrix-LHS models.
