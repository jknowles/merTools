# fastdisp: faster display of model summaries

Display model fit summary of x or x like objects, fast

## Usage

``` r
fastdisp(x, ...)

# S3 method for class 'merMod'
fastdisp(x, ...)

# S3 method for class 'merModList'
fastdisp(x, ...)
```

## Arguments

- x:

  a model object

- ...:

  additional arguments to pass to
  `arm::`[`display`](https://rdrr.io/pkg/arm/man/display.html) including
  number of digits

## Value

A list with model summary components (`call`, `coef`, `se`, `ngrps`,
`AIC`, `n`, and fit statistics), returned invisibly. The summary is also
printed to the console.

## Details

Faster than the implementation in the arm package because it avoids
refitting

The time saving is only noticeable for large, time-consuming (g)lmer
fits.

## See also

[`display`](https://rdrr.io/pkg/arm/man/display.html)

## Examples

``` r
# \donttest{
#Compare the time for displaying this modest model
require(arm)
#> Loading required package: arm
#> Loading required package: MASS
#> 
#> arm (Version 1.15-3, built: 2026-4-15)
#> Working directory is /home/runner/work/merTools/merTools/docs/reference
m1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
system.time(display(m1))
#> lmer(formula = y ~ lectage + studage + (1 | d) + (1 | s), data = InstEval)
#>             coef.est coef.se
#> (Intercept)  3.20     0.02  
#> lectage.L   -0.20     0.02  
#> lectage.Q    0.02     0.01  
#> lectage.C   -0.03     0.01  
#> lectage^4   -0.02     0.01  
#> lectage^5   -0.04     0.02  
#> studage.L    0.10     0.02  
#> studage.Q    0.01     0.02  
#> studage.C    0.02     0.02  
#> 
#> Error terms:
#>  Groups   Name        Std.Dev.
#>  s        (Intercept) 0.33    
#>  d        (Intercept) 0.52    
#>  Residual             1.18    
#> ---
#> number of obs: 73421, groups: s, 2972; d, 1128
#> AIC = 237675, DIC = 237532.5
#> deviance = 237591.5 
#>    user  system elapsed 
#>   3.399   4.307   2.020 
system.time(fastdisp(m1))
#> lmer(formula = y ~ lectage + studage + (1 | d) + (1 | s), data = InstEval)
#>             coef.est coef.se
#> (Intercept)  3.20     0.02  
#> lectage.L   -0.20     0.02  
#> lectage.Q    0.02     0.01  
#> lectage.C   -0.03     0.01  
#> lectage^4   -0.02     0.01  
#> lectage^5   -0.04     0.02  
#> studage.L    0.10     0.02  
#> studage.Q    0.01     0.02  
#> studage.C    0.02     0.02  
#> 
#> Error terms:
#>  Groups   Name        Std.Dev.
#>  s        (Intercept) 0.33    
#>  d        (Intercept) 0.52    
#>  Residual             1.18    
#> ---
#> number of obs: 73421, groups: s, 2972; d, 1128
#> AIC = 237675
#>    user  system elapsed 
#>   0.003   0.001   0.004 
# }
```
