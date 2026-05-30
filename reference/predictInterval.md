# Predict from merMod objects with a prediction interval

This function provides a way to capture model uncertainty in predictions
from multi-level models fit with `lme4`. By drawing a sampling
distribution for the random and the fixed effects and then estimating
the fitted value across that distribution, it is possible to generate a
prediction interval for fitted values that includes all variation in the
model except for variation in the covariance parameters, theta. This is
a much faster alternative than bootstrapping for models fit to medium to
large datasets.

## Usage

``` r
predictInterval(
  merMod,
  newdata = NULL,
  which = c("full", "fixed", "random", "all"),
  level = 0.8,
  n.sims = 1000,
  stat = c("median", "mean"),
  type = c("linear.prediction", "probability"),
  include.resid.var = TRUE,
  returnSims = FALSE,
  seed = NULL,
  .parallel = FALSE,
  fix.intercept.variance = FALSE,
  ignore.fixed.terms = NULL,
  new.levels = c("zero", "draw")
)
```

## Arguments

- merMod:

  a merMod object from lme4

- newdata:

  a data.frame of new data to predict

- which:

  a character specifying what to return, by default it returns the full
  interval, but you can also select to return only the fixed variation
  or the random component variation. If full is selected the resulting
  data.frame will be `nrow(newdata) * number of model levels` long

- level:

  the width of the prediction interval

- n.sims:

  number of simulation samples to construct

- stat:

  take the median or mean of simulated intervals

- type:

  type of prediction to develop

- include.resid.var:

  logical, include or exclude the residual variance for linear models

- returnSims:

  logical, should all n.sims simulations be returned?

- seed:

  numeric, optional argument to set seed for simulations

- .parallel:

  logical, use parallel computation (default FALSE)

- fix.intercept.variance:

  logical; should the variance of the intercept term be adjusted
  downwards to roughly correct for its covariance with the random
  effects, as if all the random effects are intercept effects?

- ignore.fixed.terms:

  a numeric or string vector of indexes or names of fixed effects which
  should be considered as fully known (zero variance).

- new.levels:

  character, how to treat grouping levels in `newdata` that were not
  present when the model was fit. `"zero"` (the default and the
  historical behavior) drops the random effect for such levels, so the
  prediction rests on the fixed effects plus residual variation.
  `"draw"` instead samples each unobserved group's effect from the
  estimated random-effect covariance (`VarCorr`), so the interval
  reflects between-group uncertainty – the analogue of
  `brms::posterior_predict(allow_new_levels = TRUE)`. Observations that
  share an unobserved level share the same sampled effect.

## Value

a data.frame with three columns: fit, lwr and upr. If \`returnSims\` is
TRUE the attribute \`sim.results\` contains the full simulation array.
