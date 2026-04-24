# Summarise prediction simulations into intervals

Takes the simulated prediction array and computes summary statistics
(fit, lower bound, upper bound) based on the specified statistic type
and confidence level. Handles probability transformation and the \`which
= "all"\` case with component-wise summaries.

## Usage

``` r
summarise_predictions(
  yhat_arr,
  level,
  stat.type,
  predict.type,
  N,
  merMod = NULL,
  which.eff = "full",
  pi.comps = NULL,
  is.glmm.with.response = FALSE
)
```

## Arguments

- yhat_arr:

  matrix of simulated predictions with dimensions (N x n.sims), or for
  \`which = "all"\`, the combined yhat from \`combine_components()\`.

- level:

  numeric confidence level (e.g., 0.8 for 80

- stat.type:

  character either "median" or "mean" for the central estimate.

- predict.type:

  character either "linear.prediction" or "probability".

- N:

  integer number of observations in the original newdata.

- merMod:

  merMod object, required when \`predict.type = "probability"\` to
  access the link function.

- which.eff:

  character the effect type: "full", "fixed", "random", or "all".

- pi.comps:

  list of component prediction matrices (from \`combine_components()\`
  when \`which = "all"\`), or NULL for other cases.

## Value

data.frame with columns \`fit\`, \`upr\`, \`lwr\`. For \`which =
"all"\`, includes additional columns \`effect\` and \`obs\`.
