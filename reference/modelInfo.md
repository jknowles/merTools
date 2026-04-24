# Extract model information from a merMod

Extract model information from a merMod

## Usage

``` r
modelInfo(object)
```

## Arguments

- object:

  a merMod object

## Value

Simple summary information about the object, number of observations,
number of grouping terms, AIC, and residual standard deviation

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelInfo(mod[[1]])
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
lapply(mod, modelInfo)
#> [[1]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[2]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[3]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[4]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[5]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[6]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[7]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[8]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[9]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
#> [[10]]
#>   n.obs n.lvls      AIC   sigma
#> 1   180      1 1755.628 25.5918
#> 
# }
```
