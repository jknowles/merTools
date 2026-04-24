# Extract theta parameters from a merMod model

A convenience function that returns the theta parameters for a
[`merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html) object.

## Usage

``` r
thetaExtract(merMod)
```

## Arguments

- merMod:

  a valid merMod object

## Value

a vector of the covariance, theta, parameters from a
[`merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html)

## See also

merMod

## Examples

``` r
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
thetaExtract(fm1) #(a numeric vector of the covariance parameters)
#> [1] 0.96674177 0.01516906 0.23090995
```
