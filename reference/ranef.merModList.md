# Extract random-effects estimates for a merModList

Extract random-effects estimates for a merModList

## Usage

``` r
# S3 method for class 'merModList'
ranef(object, ...)
```

## Arguments

- object:

  a [`merModList`](merModList-class.md) list fitted models with random
  effects

- ...:

  some methods for these generic functions require additional arguments.

## Value

a named, numeric vector of random-effects estimates.

## Details

Extract the estimates of the random-effects parameters from a list of
fitted [`merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html)
models. Takes the mean of the individual `ranef` objects for each of the
component models in the [`merModList`](merModList-class.md).

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
ranef(mod)
#> $Subject
#>     (Intercept)        Days
#> 308  0.22585510  0.91989758
#> 309 -4.03987381 -0.86196806
#> 310 -3.89604090 -0.54488565
#> 330  2.36906196 -0.48143503
#> 331  2.22603126 -0.30699116
#> 332  0.90395679 -0.02721770
#> 333  1.68405086 -0.02236361
#> 334 -0.72326151  0.10745816
#> 335 -0.03336684 -1.07521652
#> 337  3.48904868  0.86282652
#> 349 -2.52102286  0.11734322
#> 350 -1.30700342  0.66142178
#> 351  0.45778642 -0.30152621
#> 352  2.08636782  0.35360011
#> 369  0.32754656  0.08722149
#> 370 -2.56129993  0.48224850
#> 371  0.08070461 -0.09881562
#> 372  1.23145921  0.12840221
#> 
# }
```
