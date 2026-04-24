# Extract data.frame of random effect statistics from merMod List

Extract data.frame of random effect statistics from merMod List

## Usage

``` r
modelRandEffStats(modList)
```

## Arguments

- modList:

  a list of multilevel models

## Value

a data.frame

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelRandEffStats(mod)
#>                    term    group    estimate    std.error
#> 1 cor__(Intercept).Days  Subject  0.06555124 2.733082e-11
#> 2       sd__(Intercept)  Subject 24.74065799 8.339641e-10
#> 3              sd__Days  Subject  5.92213766 6.571468e-11
#> 4       sd__Observation Residual 25.59179572 8.821610e-11
# }
```
