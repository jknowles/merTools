# Extracts random effects

Extracts random effect terms from an lme4 model

## Usage

``` r
REextract(merMod)
```

## Arguments

- merMod:

  a merMod object from the lme4 package

## Value

a data frame with the following columns

- groupFctr:

  The name of the grouping factor associated with the random effects

- groupID:

  The level of the grouping factor associated with the random effects

- 'term':

  One column per random effect, the name is derived from the merMod

- 'term'\_se:

  One column per random effect, the name is derived from the merMod

## Examples

``` r
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
rfx <- REextract(m2)
#Note the column names
head(rfx)
#>     groupFctr groupID (Intercept)      Days (Intercept)_se  Days_se
#> 308   Subject     308    2.258551  9.198976       12.07086 2.304839
#> 309   Subject     309  -40.398738 -8.619681       12.07086 2.304839
#> 310   Subject     310  -38.960409 -5.448856       12.07086 2.304839
#> 330   Subject     330   23.690620 -4.814350       12.07086 2.304839
#> 331   Subject     331   22.260313 -3.069912       12.07086 2.304839
#> 332   Subject     332    9.039568 -0.272177       12.07086 2.304839
```
