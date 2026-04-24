# Extract the standard deviation of the random effects from a merMod object

Extract the standard deviation of the random effects from a merMod
object

## Usage

``` r
REsdExtract(model)
```

## Arguments

- model:

  an object that inherits from class merMod

## Value

a numeric vector for standard deviations of the random effects

## Examples

``` r
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REsdExtract(fm1)
#> Subject.(Intercept)        Subject.Days 
#>           24.740658            5.922138 
```
