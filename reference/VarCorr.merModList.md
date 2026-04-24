# Extract the variances and correlations for random effects from a merMod list

Extract the variances and correlations for random effects from a merMod
list

## Usage

``` r
# S3 method for class 'merModList'
VarCorr(x, sigma = 1, rdig = 3L, ...)
```

## Arguments

- x:

  a [`merModList`](merModList-class.md) list fitted models with random
  effects

- sigma:

  an optional numeric value used as a multiplier for the standard
  deviations.

- rdig:

  the number of digits to round to, integer

- ...:

  Ignored for the `as.data.frame` method; passed to other
  [`print()`](https://rdrr.io/r/base/print.html) methods for the
  [`print()`](https://rdrr.io/r/base/print.html) method.

## Value

a list with two elements "stddev" and "correlation" for the standard
deviations and correlations averaged across models in the list

## Examples

``` r
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
VarCorr(mod)
#> $stddev
#> $stddev$Subject
#> (Intercept)        Days 
#>   24.740658    5.922138 
#> 
#> 
#> $correlation
#> $correlation$Subject
#>             (Intercept)       Days
#> (Intercept)  1.00000000 0.06555124
#> Days         0.06555124 1.00000000
#> 
#> 
```
