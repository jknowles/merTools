# Extract fixed-effects estimates for a merModList

Extract fixed-effects estimates for a merModList

## Usage

``` r
# S3 method for class 'merModList'
fixef(object, add.dropped = FALSE, ...)
```

## Arguments

- object:

  any fitted model object from which fixed effects estimates can be
  extracted.

- add.dropped:

  for models with rank-deficient design matrix, reconstitute the
  full-length parameter vector by adding `NA` values in appropriate
  locations?

- ...:

  optional additional arguments. Currently none are used in any methods.

## Value

a named, numeric vector of fixed-effects estimates.

## Details

Extract the estimates of the fixed-effects parameters from a list of
fitted `merMod` models. Takes the mean of the individual `fixef` objects
for each of the component models in the `merModList`.

## Examples

``` r
# \donttest{
sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
fixef(mod)
#> (Intercept)        Days 
#>   251.40510    10.46729 
# }
```
