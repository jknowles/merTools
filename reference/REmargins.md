# Calculate the predicted value for each observation across the distribution of the random effect terms.

`REmargins` calculates the average predicted value for each row of a new
data frame across the distribution of
[`expectedRank`](https://jknowles.github.io/merTools/reference/expectedRank.md)
for a merMod object. This allows the user to make meaningful comparisons
about the influence of random effect terms on the scale of the response
variable, for user-defined inputs, and accounting for the variability in
grouping terms.

## Usage

``` r
REmargins(
  merMod,
  newdata = NULL,
  groupFctr = NULL,
  term = NULL,
  breaks = 4,
  .parallel = FALSE,
  ...
)
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
  by `ranef`. If the length is \> 1 then the combined effect of all
  listed groups will calculated and marginalized over co-occurrences of
  those groups if desired.

- term:

  The name of the random coefficient of interest. This is the variable
  to the left of the pipe, `|`, in the \[g\]lmer formula. Partial
  matching is attempted on the intercept term so the following character
  strings will all return rankings based on the intercept (*provided
  that they do not match the name of another random coefficient for that
  factor*): `c("(Intercept)", "Int", "intercep", ...)`.

- breaks:

  an integer representing the number of bins to divide the group effects
  into, the default is 3.

- .parallel, :

  logical should parallel computation be used, default is TRUE

- ...:

  additional arguments to pass to
  [`predictInterval`](https://jknowles.github.io/merTools/reference/predictInterval.md)

## Value

A data.frame with all unique combinations of the number of cases, rows
in the newdata element:

- ...:

  The columns of the original data taken from `newdata`

- case:

  The row number of the observation from newdata. Each row in newdata
  will be repeated for all unique levels of the grouping_var, term, and
  breaks.

- grouping_var:

  The grouping variable the random effect is being marginalized over.

- term:

  The term for the grouping variable the random effect is being
  marginalized over.

- breaks:

  The ntile of the effect size for `grouping_var` and `term`

- original_group_level:

  The original grouping value for this `case`

- fit_combined:

  The predicted value from `predictInterval` for this case simulated at
  the Nth ntile of the expected rank distribution of `grouping_var` and
  `term`

- upr_combined:

  The upper bound of the predicted value.

- lwr_combined:

  The lower bound of the predicted value.

- fit_XX:

  For each grouping term in newdata the predicted value is decomposed
  into its fit components via predictInterval and these are all returned
  here

- upr_XX:

  The upper bound for the effect of each grouping term

- lwr_XX:

  The lower bound for the effect of each grouping term

- fit_fixed:

  The predicted fit with all the grouping terms set to 0 (average)

- upr_fixed:

  The upper bound fit with all the grouping terms set to 0 (average)

- lwr_fixed:

  The lower bound fit with all the grouping terms set to 0 (average)

## Details

The function simulates the

The function predicts the response at every level in the random effect
term specified by the user. Then, the expected rank of each group level
is binned to the number of bins specified by the user. Finally, a
weighted mean of the fitted value for all observations in each bin of
the expected ranks is calculated using the inverse of the variance as
the weight – so that less precise estimates are downweighted in the
calculation of the mean for the bin. Finally, a standard error for the
bin mean is calculated.

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
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:10,])

# You can also pass additional arguments to predictInterval through REimpact
 g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
 margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("s"),
                        breaks = 4)
 margin_df <- REmargins(g1, newdata = InstEval[20:25, ], groupFctr = c("d"),
                         breaks = 3)
# }
```
