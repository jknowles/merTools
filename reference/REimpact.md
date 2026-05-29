# Calculate the weighted mean of fitted values for various levels of random effect terms.

`REimpact` calculates the average predicted value for each row of a new
data frame across the distribution of
[`expectedRank`](https://jknowles.github.io/merTools/reference/expectedRank.md)
for a merMod object. This allows the user to make meaningful comparisons
about the influence of random effect terms on the scale of the response
variable, for user-defined inputs, and accounting for the variability in
grouping terms.

## Usage

``` r
REimpact(merMod, newdata, groupFctr = NULL, term = NULL, breaks = 3, ...)
```

## Arguments

- merMod:

  An object of class merMod

- newdata:

  a data frame of observations to calculate group-level differences for

- groupFctr:

  The name of the grouping factor over which the random coefficient of
  interest varies. This is the variable to the right of the pipe, `|`,
  in the \[g\]lmer formula. This parameter is optional, if not
  specified, it will perform the calculation for the first effect listed
  by `ranef`.

- term:

  The name of the random coefficient of interest. This is the variable
  to the left of the pipe, `|`, in the \[g\]lmer formula. Partial
  matching is attempted on the intercept term so the following character
  strings will all return rankings based on the intercept (*provided
  that they do not match the name of another random coefficient for that
  factor*): `c("(Intercept)", "Int", "intercep", ...)`.

- breaks:

  an integer representing the number of bins to divide the group effects
  into, the default is 3; alternatively it can specify breaks from 0-100
  for how to cut the expected rank distribution

- ...:

  additional arguments to pass to
  [`predictInterval`](https://jknowles.github.io/merTools/reference/predictInterval.md)

## Value

A data.frame with all unique combinations of the number of cases, rows
in the newdata element, and number of bins:

- case:

  The row number of the observation from newdata.

- bin:

  The ranking bin for the expected rank, the higher the bin number, the
  greater the expected rank of the groups in that bin.

- AvgFitWght:

  The weighted mean of the fitted values for case i in bin k

- AvgFitWghtSE:

  The standard deviation of the mean of the fitted values for case i in
  bin k.

- nobs:

  The number of group effects contained in that bin.

## Details

The function predicts the response at every level in the random effect
term specified by the user. Then, the expected rank of each group level
is binned to the number of bins specified by the user. Finally, a
weighted mean of the fitted value for all observations in each bin of
the expected ranks is calculated using the inverse of the variance as
the weight – so that less precise estimates are downweighted in the
calculation of the mean for the bin. Finally, a standard error for the
bin mean is calculated.

This function uses the formula for variance of a weighted mean
recommended by Cochran (1977).

## References

Gatz, DF and Smith, L. The Standard Error of a Weighted Mean
Concentration. I. Bootstrapping vs other methods. *Atmospheric
Environment*. 1995;11(2)1185-1193. Available at
<https://www.sciencedirect.com/science/article/pii/135223109400210C>

Cochran, WG. 1977. Sampling Techniques (3rd Edition). Wiley, New York.

## See also

[`expectedRank`](https://jknowles.github.io/merTools/reference/expectedRank.md),
[`predictInterval`](https://jknowles.github.io/merTools/reference/predictInterval.md)

## Examples

``` r
# \donttest{
#For a one-level random intercept model
m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
m1.er <- REimpact(m1, newdata = sleepstudy[1, ], breaks = 2)
#For a one-level random intercept model with multiple random terms
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#ranked by the random slope on Days
m2.er1 <- REimpact(m2,  newdata = sleepstudy[1, ],
           groupFctr = "Subject", term="Days")
#ranked by the random intercept
m2.er2 <- REimpact(m2, newdata = sleepstudy[1, ],
             groupFctr = "Subject", term="int")

# You can also pass additional arguments to predictInterval through REimpact
g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
zed <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", n.sims = 50,
                include.resid.var = TRUE)
#> Warning: executing %dopar% sequentially: no parallel backend registered
zed2 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "s", n.sims = 50,
                 include.resid.var = TRUE)
zed3 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", breaks = 5,
                n.sims = 50, include.resid.var = TRUE)
# }
```
