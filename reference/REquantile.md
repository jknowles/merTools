# Identify group level associated with RE quantile

For a user specified quantile (or quantiles) of the random effect terms
in a merMod object. This allows the user to easily identify the
observation associated with the nth percentile effect.

## Usage

``` r
REquantile(merMod, quantile, groupFctr, term = "(Intercept)")
```

## Arguments

- merMod:

  a merMod object with one or more random effect levels

- quantile:

  a numeric vector with values between 0 and 100 for quantiles

- groupFctr:

  a character of the name of the random effect grouping factor to
  extract quantiles from

- term:

  a character of the random effect to extract for the grouping factor
  specified. Default is the intercept.

## Value

a vector of the level of the random effect grouping term that
corresponds to each quantile

## Examples

``` r
# \donttest{
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REquantile(fm1, quantile = 0.25, groupFctr = "Subject")
#> Number of observations < 20, random effect quantiles may not be well-defined.
#> [1] "349"
REquantile(fm1, quantile = 0.25, groupFctr = "Subject", term = "Days")
#> Number of observations < 20, random effect quantiles may not be well-defined.
#> [1] "330"
# }
```
