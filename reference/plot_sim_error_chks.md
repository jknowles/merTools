# Extract all warning msgs from a merMod object

Extract all warning msgs from a merMod object

## Usage

``` r
plot_sim_error_chks(
  type = c("FE", "RE"),
  level = 0.95,
  stat = c("mean", "median"),
  sd = TRUE,
  sigmaScale = NULL,
  oddsRatio = FALSE,
  labs = FALSE,
  facet = TRUE
)
```

## Arguments

- type:

  check a fixed or random effect

- level:

  the width of the confidence interval

- stat:

  a character value indicating the variable name in data of the midpoint
  of the estimated interval, e.g. "mean" or "median"

- sd:

  a logical indicating whether or not to plot error bars around the
  estimates (default is TRUE). Calculates the width of the error bars
  based on `level` and the variable named "sd" in `data`

- sigmaScale:

  a numeric value to divide the estimate and the standard deviation by
  in the case of doing an effect size calculation

- oddsRatio:

  logical, should the parameters be converted to odds ratios before
  plotting

- labs:

  logical, include the labels of the groups on the x-axis

- facet:

  Accepts either logical (`TRUE`) or `list` to specify which random
  effects to plot. If `TRUE`, facets by both `groupFctr` and `term`. If
  `list` selects the panel specified by the named elements of the list
