# Bootstrap a subset of an lme4 model

Bootstrap a subset of an lme4 model

## Usage

``` r
subBoot(merMod, n = NULL, FUN, R = 100, seed = NULL, warn = FALSE)
```

## Arguments

- merMod:

  a valid merMod object

- n:

  the number of rows to sample from the original data in the merMod
  object, by default will resample the entire model frame

- FUN:

  the function to apply to each bootstrapped model

- R:

  the number of bootstrap replicates, default is 100

- seed:

  numeric, optional argument to set seed for simulations

- warn:

  logical, if TRUE, warnings from lmer will be issued, otherwise they
  will be suppressed default is FALSE

## Value

a data.frame of parameters extracted from each of the R replications.
The original values are appended to the top of the matrix.

## Details

This function allows users to estimate parameters of a large merMod
object using bootstraps on a subset of the data.

## Examples

``` r
# \donttest{
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: Reaction ~ Days + (Days | Subject)
#>    Data: sleepstudy
#> REML criterion at convergence: 1743.628
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  Subject  (Intercept) 24.741        
#>           Days         5.922   0.07 
#>  Residual             25.592        
#> Number of obs: 180, groups:  Subject, 18
#> Fixed Effects:
#> (Intercept)         Days  
#>      251.41        10.47  
resultMatrix <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20)
#> Warnings set to off by default, not all submodels may have converged.
# }
```
