# Extract averaged fixed effect parameters across a list of merMod objects

Extract averaged fixed effect parameters across a list of merMod objects

## Usage

``` r
modelFixedEff(modList, ...)
```

## Arguments

- modList:

  an object of class merModList

- ...:

  additional arguments to pass to
  [`tidy`](https://generics.r-lib.org/reference/tidy.html)

## Value

a data.frame of the averaged fixed effect parameters

## Details

The Rubin correction for combining estimates and standard errors from
Rubin (1987) is applied to adjust for the within and between imputation
variances.

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelFixedEff(mod)
#> Warning: Between imputation variance is very small, are imputation sets too similar?
#>          term  estimate std.error statistic           df
#> 1 (Intercept) 251.40510  6.824597 36.838090 6.956236e+49
#> 2        Days  10.46729  1.545790  6.771481 4.537864e+51
# }
```
