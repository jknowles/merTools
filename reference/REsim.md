# Simulate random effects from merMod `REsim` simulates random effects from merMod object posterior distributions

Simulate random effects from merMod `REsim` simulates random effects
from merMod object posterior distributions

## Usage

``` r
REsim(merMod, n.sims = 200, oddsRatio = FALSE, seed = NULL)
```

## Arguments

- merMod:

  a merMod object from the lme4 package

- n.sims:

  number of simulations to use

- oddsRatio:

  logical, should parameters be converted to odds ratios?

- seed:

  numeric, optional argument to set seed for simulations

## Value

a data frame with the following columns

- `groupFctr`:

  Name of the grouping factor

- `groupID`:

  Level of the grouping factor

- `term`:

  Name of random term (intercept/coefficient)

- `mean`:

  Mean of the simulations

- `median`:

  Median of the simulations

- `sd`:

  Standard deviation of the simulations, `NA` if `oddsRatio=TRUE`

## Details

Use the Gelman sim technique to build empirical Bayes estimates. Uses
the sim function in the arm package

## Examples

``` r
require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
re2 <- REsim(m2, 25)
head(re2)
#>   groupFctr groupID        term       mean     median       sd
#> 1   Subject     308 (Intercept)   5.618048   4.175001 10.53275
#> 2   Subject     309 (Intercept) -39.389741 -38.657014 15.41677
#> 3   Subject     310 (Intercept) -34.562139 -37.224606 12.07460
#> 4   Subject     330 (Intercept)  20.392765  20.362758 12.55282
#> 5   Subject     331 (Intercept)  24.154543  22.841477 11.93401
#> 6   Subject     332 (Intercept)   6.339342   4.354799 11.50272
```
