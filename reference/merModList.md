# Apply a multilevel model to a list of data frames

Apply a multilevel model to a list of data frames

Apply a Bayesian multilevel model to a list of data frames

Apply a generalized linear multilevel model to a list of data frames

Apply a Bayesian generalized linear multilevel model to a list of data
frames

## Usage

``` r
lmerModList(formula, data, parallel = FALSE, ...)

blmerModList(formula, data, parallel = FALSE, ...)

glmerModList(formula, data, parallel = FALSE, ...)

bglmerModList(formula, data, parallel = FALSE, ...)
```

## Arguments

- formula:

  a formula to pass through compatible with merMod

- data:

  a list object with each element being a data.frame

- parallel:

  logical, should the models be run in parallel? Default FALSE. If so,
  the \`future_lapply\` function from the \`future.apply\` package is
  used. See details.

- ...:

  additional arguments to pass to the estimating function

## Value

a list of fitted merMod objects of class merModList

a merModList

a merModList

a merModList

## Details

Parallel computing is provided by the \`futures\` package, and its
extension the \`future.apply\` package to provide the \`future_lapply\`
function for easy parallel computations on lists. To use this package,
simply register a parallel backend using the \`plan()\` function from
\`futures\` - an example is to use \`plan(multisession)\`

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
#> (Intercept)  251.405     6.825    36.838 1.312371e+50
#> Days          10.467     1.546     6.771 2.257823e+51
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
