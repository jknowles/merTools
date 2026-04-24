# Summarize a merMod list

Summarize a merMod list

## Usage

``` r
# S3 method for class 'merModList'
print(x, ...)
```

## Arguments

- x:

  a modList of class merModList

- ...:

  additional arguments

## Value

a summary object of model information

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
summary(mod)
#> Warning: Between imputation variance is very small, are imputation sets too similar?
#> [1] "Linear mixed model fit by REML"
#> Model family: 
#> lmer(formula = Reaction ~ Days + (Days | Subject), data = d)
#> 
#> Fixed Effects:
#>             estimate std.error statistic           df
#> (Intercept)  251.405     6.825    36.838 1.018856e+50
#> Days          10.467     1.546     6.771 3.794005e+51
#> 
#> Random Effects:
#> 
#> Error Term Standard Deviations by Level:
#> 
#> Subject
#> (Intercept)        Days 
#>      24.741       5.922 
#> 
#> 
#> Error Term Correlations:
#> 
#> Subject
#>             (Intercept) Days 
#> (Intercept) 1.000       0.066
#> Days        0.066       1.000
#> 
#> 
#> Residual Error = 25.592 
#> 
#> ---Groups
#> number of obs: 180, groups: Subject, 18
#> 
#> Model Fit Stats
#> AIC = 1755.6
#> Residual standard deviation = 25.592 
# }
```
