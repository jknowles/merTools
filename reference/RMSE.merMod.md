# Estimate the Root Mean Squared Error (RMSE) for a lmerMod

Extract the Root Mean Squared Error for a lmerMod object

## Usage

``` r
RMSE.merMod(merMod, scale = FALSE)
```

## Arguments

- merMod:

  a lmerMod object from the lme4 package

- scale:

  logical, should the result be returned on the scale of response
  variable standard deviations?

## Value

a numeric which represents the RMSE

## Examples

``` r
require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
RMSE.merMod(m2)
#> [1] 23.43805
```
