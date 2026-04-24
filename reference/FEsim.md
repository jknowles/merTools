# Simulate fixed effects from merMod `FEsim` simulates fixed effects from merMod object posterior distributions

Simulate fixed effects from merMod `FEsim` simulates fixed effects from
merMod object posterior distributions

## Usage

``` r
FEsim(merMod, n.sims = 200, oddsRatio = FALSE, seed = NULL)
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

- `term`:

  Name of fixed term (intercept/coefficient)

- `mean`:

  Mean of the simulations

- `median`:

  Median of the simulations

- `sd`:

  Standard deviation of the simulations, `NA` if `oddsRatio=TRUE`

## Details

Use the Gelman sim technique to build fixed effect estimates and
confidence intervals. Uses the sim function in the arm package

## Examples

``` r
require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fe2 <- FEsim(m2, 25)
head(fe2)
#>          term      mean    median       sd
#> 1 (Intercept) 250.17301 250.48344 5.332610
#> 2        Days  10.60686  10.71097 1.305392
```
