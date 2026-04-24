# Extract the correlations between the slopes and the intercepts from a model

Extract the correlations between the slopes and the intercepts from a
model

## Usage

``` r
REcorrExtract(model)
```

## Arguments

- model:

  an object that inherits from class merMod

## Value

a numeric vector of the correlations among the effects

## Examples

``` r
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REcorrExtract(fm1)
#> [1] 0.06555124
```
