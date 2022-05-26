[![Travis-CI Build
Status](https://travis-ci.org/jknowles/merTools.png?branch=master)](https://travis-ci.org/jknowles/merTools)
[![Coverage
Status](https://coveralls.io/repos/jknowles/merTools/badge.svg?branch=master)](https://coveralls.io/r/jknowles/merTools?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/merTools)](https://cran.r-project.org/package=merTools)
[![Downloads](http://cranlogs.r-pkg.org/badges/merTools)](https://cran.r-project.org/package=merTools)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/merTools)](https://cran.r-project.org/package=merTools)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# merTools

A package for getting the most out of large multilevel models in R

by Jared E. Knowles and Carl Frederick

Working with generalized linear mixed models (GLMM) and linear mixed
models (LMM) has become increasingly easy with advances in the `lme4`
package. As we have found ourselves using these models more and more
within our work, we, the authors, have developed a set of tools for
simplifying and speeding up common tasks for interacting with `merMod`
objects from `lme4`. This package provides those tools.

## Installation

``` r
# development version
library(devtools)
install_github("jknowles/merTools")

# CRAN version
install.packages("merTools")
```

## Recent Updates

### merTools 0.5.0

#### New Features

-   `subBoot` now works with `glmerMod` objects as well
-   `reMargins` a new function that allows the user to marginalize the
    prediction over breaks in the distribution of random effect
    distributions, see `?reMargins` and the new `reMargins` vignette
    (closes \#73)

#### Bug fixes

-   Fixed an issue where known convergence errors were issuing warnings
    and causing the test suite to not work
-   Fixed an issue where models with a random slope, no intercept, and
    no fixed term were unable to be predicted (#101)
-   Fixed an issue with shinyMer not working with substantive fixed
    effects (#93)

### merTools 0.4.1

#### New Features

-   Standard errors reported by `merModList` functions now apply the
    Rubin correction for multiple imputation

#### Bug fixes

-   Contribution by Alex Whitworth (@alexWhitworth) adding error
    checking to plotting functions

### merTools 0.4.0

#### New Features

-   Added vignette on using multilevel models with multiply imputed data
-   Added `fixef` and `ranef` generics for `merModList` objects
-   Added `fastdisp` generic for `merModList`
-   Added `summary` generic for `merModList`
-   Added `print` generic for `merModList`
-   Documented all generics for `merModList` including examples and a
    new imputation vignette
-   Added `modelInfo` generic for `merMod` objects that provides simple
    summary stats about a whole model

#### Bug Fixes

-   Fix bug that returned NaN for `std.error` of a multiply imputed
    `merModList` when calling `modelRandEffStats`
-   Fixed bug in `REimpact` where some column names in `newdata` would
    prevent the prediction intervals from being computed correctly.
    Users will now be warned.
-   Fixed bug in `wiggle` where documentation incorrectly stated the
    arguments to the function and the documentation did not describe
    function correctly

See [NEWS.md](https://github.com/jknowles/merTools/blob/master/NEWS.md)
for more details.

## Shiny App and Demo

The easiest way to demo the features of this application is to use the
bundled Shiny application which launches a number of the metrics here to
aide in exploring the model. To do this:

``` r
library(merTools)
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data
```

![](man/figures/README-predPanel.png)

On the first tab, the function presents the prediction intervals for the
data selected by user which are calculated using the `predictInterval`
function within the package. This function calculates prediction
intervals quickly by sampling from the simulated distribution of the
fixed effect and random effect terms and combining these simulated
estimates to produce a distribution of predictions for each observation.
This allows prediction intervals to be generated from very large models
where the use of `bootMer` would not be feasible computationally.

![](man/figures/README-effPanel.png)

On the next tab the distribution of the fixed effect and group-level
effects is depicted on confidence interval plots. These are useful for
diagnostics and provide a way to inspect the relative magnitudes of
various parameters. This tab makes use of four related functions in
`merTools`: `FEsim`, `plotFEsim`, `REsim` and `plotREsim` which are
available to be used on their own as well.

![](man/figures/README-substPanel.png)

On the third tab are some convenient ways to show the influence or
magnitude of effects by leveraging the power of `predictInterval`. For
each case, up to 12, in the selected data type, the user can view the
impact of changing either one of the fixed effect or one of the grouping
level terms. Using the `REimpact` function, each case is simulated with
the model’s prediction if all else was held equal, but the observation
was moved through the distribution of the fixed effect or the random
effect term. This is plotted on the scale of the dependent variable,
which allows the user to compare the magnitude of effects across
variables, and also between models on the same data.

## Predicting

Standard prediction looks like so.

``` r
predict(m1, newdata = InstEval[1:10, ])
#>        1        2        3        4        5        6        7        8 
#> 3.146337 3.165212 3.398499 3.114249 3.320686 3.252670 4.180897 3.845219 
#>        9       10 
#> 3.779337 3.331013
```

With `predictInterval` we obtain predictions that are more like the
standard objects produced by `lm` and `glm`:

``` r
#predictInterval(m1, newdata = InstEval[1:10, ]) # all other parameters are optional
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 500, level = 0.9, 
                stat = 'median')
#>         fit      upr      lwr
#> 1  3.132412 5.049890 1.150327
#> 2  3.069004 4.967174 1.349411
#> 3  3.371890 5.225487 1.203104
#> 4  3.065470 5.098975 1.052925
#> 5  3.294912 5.262798 1.162989
#> 6  3.244553 5.288931 1.138605
#> 7  4.231166 6.258813 2.424893
#> 8  3.781273 5.604849 2.099314
#> 9  3.836202 5.747436 1.676611
#> 10 3.367951 5.276930 1.423142
```

Note that `predictInterval` is slower because it is computing
simulations. It can also return all of the simulated `yhat` values as an
attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily
to draw the distributions of the parameters of the model. It then
combines these simulated values to create a distribution of the `yhat`
for each observation.

### Inspecting the Prediction Components

We can also explore the components of the prediction interval by asking
`predictInterval` to return specific components of the prediction
interval.

``` r
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 200, level = 0.9, 
                stat = 'median', which = "all")
#>      effect           fit      upr        lwr obs
#> 1  combined  3.2540525131 4.851934  1.2475328   1
#> 2  combined  3.1753216228 5.224952  1.1905459   2
#> 3  combined  3.4454221587 5.532090  1.2890110   3
#> 4  combined  3.1578660116 5.317547  1.1327241   4
#> 5  combined  3.4784433122 5.479333  1.4475864   5
#> 6  combined  3.1395561625 5.277999  1.2131769   6
#> 7  combined  4.0953033332 5.897310  2.2026540   7
#> 8  combined  3.8497708965 5.520022  1.8972317   8
#> 9  combined  3.8399908781 5.324169  1.6335418   9
#> 10 combined  3.3522919893 5.138446  1.7777244  10
#> 11        s  0.1609574017 2.280438 -1.6557246   1
#> 12        s  0.3452691401 2.398365 -1.6241005   2
#> 13        s -0.0006663979 2.062559 -1.6416237   3
#> 14        s  0.0210520917 1.942374 -1.6583963   4
#> 15        s -0.1742755831 2.359372 -2.0115378   5
#> 16        s  0.0392215140 2.138145 -1.9365056   6
#> 17        s  0.3473602858 2.382207 -1.5927409   7
#> 18        s  0.1321817157 2.040468 -1.5521957   8
#> 19        s  0.3581967760 2.301422 -1.6630743   9
#> 20        s  0.2418774146 2.172069 -1.6219465  10
#> 21        d -0.2188230250 1.578214 -2.0762967   1
#> 22        d -0.3078932312 1.808852 -2.0341853   2
#> 23        d  0.0392979203 1.911539 -1.8048353   3
#> 24        d -0.0656348864 1.698348 -2.2446730   4
#> 25        d  0.1754887186 2.118397 -2.0089090   5
#> 26        d  0.1150995845 2.143500 -2.0168937   6
#> 27        d  0.7218818083 2.389070 -1.0393478   7
#> 28        d  0.1930955117 2.138620 -1.9178248   8
#> 29        d  0.3260866954 2.507140 -1.8411643   9
#> 30        d -0.3521779048 1.649898 -2.3212154  10
#> 31    fixed  3.2538642529 5.066044  1.3034822   1
#> 32    fixed  3.1234760300 4.973536  0.8102403   2
#> 33    fixed  3.1557416254 5.206102  1.3762575   3
#> 34    fixed  3.1122406870 5.091127  1.0914009   4
#> 35    fixed  3.4325167242 5.221134  1.2653626   5
#> 36    fixed  3.3740892423 5.183653  1.2385467   6
#> 37    fixed  3.0971849967 5.158864  1.4192060   7
#> 38    fixed  3.1386171238 4.930877  1.2996725   8
#> 39    fixed  3.1457732415 5.370632  1.5128355   9
#> 40    fixed  3.1390823897 5.329963  1.5155335  10
```

This can lead to some useful plotting:

``` r
library(ggplot2)
plotdf <- predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 2000, 
                          level = 0.9, stat = 'median', which = "all", 
                          include.resid.var = FALSE)
plotdfb <- predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 2000, 
                          level = 0.9, stat = 'median', which = "all", 
                          include.resid.var = TRUE)

plotdf <- dplyr::bind_rows(plotdf, plotdfb, .id = "residVar")
plotdf$residVar <- ifelse(plotdf$residVar == 1, "No Model Variance", 
                          "Model Variance")

ggplot(plotdf, aes(x = obs, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  scale_x_continuous(breaks = c(1, 10)) +
  facet_grid(residVar~effect) + theme_bw()
```

![](man/figures/README_unnamed-chunk-8-1.png)<!-- -->

We can also investigate the makeup of the prediction for each
observation.

``` r
ggplot(plotdf[plotdf$obs < 6,], 
       aes(x = effect, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  facet_grid(residVar~obs) + theme_bw()
```

![](man/figures/README_unnamed-chunk-9-1.png)<!-- -->

## Plotting

`merTools` also provides functionality for inspecting `merMod` objects
visually. The easiest are getting the posterior distributions of both
fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22767582  3.22761955 0.02049772
#> 2    service1 -0.07086013 -0.07125082 0.01227076
#> 3   lectage.L -0.18444289 -0.18381933 0.01468229
#> 4   lectage.Q  0.02078291  0.02141907 0.01286262
#> 5   lectage.C -0.02501968 -0.02435823 0.01221798
#> 6   lectage^4 -0.02014908 -0.02107177 0.01212984
```

And we can also plot this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](man/figures/README_FEsimPlot-1.png)<!-- -->

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean       median        sd
#> 1         s       1 (Intercept)  0.14968443  0.126556409 0.2830427
#> 2         s       2 (Intercept) -0.02687232 -0.024849270 0.2828727
#> 3         s       3 (Intercept)  0.26895812  0.252701182 0.2732403
#> 4         s       4 (Intercept)  0.19719511  0.169396113 0.2771770
#> 5         s       5 (Intercept)  0.02378876  0.003333246 0.2867102
#> 6         s       6 (Intercept)  0.08342151  0.074260609 0.2350585
```

``` r
plotREsim(REsim(m1, n.sims = 100), stat = 'median', sd = TRUE)
```

![](man/figures/README_reSimplot-1.png)<!-- -->

Note that `plotREsim` highlights group levels that have a simulated
distribution that does not overlap 0 – these appear darker. The lighter
bars represent grouping levels that are not distinguishable from 0 in
the data.

Sometimes the random effects can be hard to interpret and not all of
them are meaningfully different from zero. To help with this `merTools`
provides the `expectedRank` function, which provides the percentile
ranks for the observed groups in the random effect distribution taking
into account both the magnitude and uncertainty of the estimated effect
for each group.

``` r
ranks <- expectedRank(m1, groupFctr = "d")
head(ranks)
#>   groupFctr groupLevel      term   estimate  std.error       ER pctER
#> 2         d          1 Intercept  0.3944919 0.08665152 835.3005    74
#> 3         d          6 Intercept -0.4428949 0.03901988 239.5363    21
#> 4         d          7 Intercept  0.6562681 0.03717200 997.3569    88
#> 5         d          8 Intercept -0.6430680 0.02210017 138.3445    12
#> 6         d         12 Intercept  0.1902940 0.04024063 702.3410    62
#> 7         d         13 Intercept  0.2497464 0.03216255 750.0174    66
```

A nice features `expectedRank` is that you can return the expected rank
for all factors simultaneously and use them:

``` r
ranks <- expectedRank(m1)
head(ranks)
#>   groupFctr groupLevel      term    estimate  std.error       ER pctER
#> 2         s          1 Intercept  0.16732800 0.08165665 1931.570    65
#> 3         s          2 Intercept -0.04409538 0.09234250 1368.160    46
#> 4         s          3 Intercept  0.30382219 0.05204082 2309.941    78
#> 5         s          4 Intercept  0.24756175 0.06641699 2151.828    72
#> 6         s          5 Intercept  0.05232329 0.08174130 1627.693    55
#> 7         s          6 Intercept  0.10191653 0.06648394 1772.548    60

ggplot(ranks, aes(x = term, y = estimate)) + 
  geom_violin(fill = "gray50") + facet_wrap(~groupFctr) +
  theme_bw()
```

![](man/figures/README_unnamed-chunk-13-1.png)<!-- -->

## Effect Simulation

It can still be difficult to interpret the results of LMM and GLMM
models, especially the relative influence of varying parameters on the
predicted outcome. This is where the `REimpact` and the `wiggle`
functions in `merTools` can be handy.

``` r
impSim <- REimpact(m1, InstEval[7, ], groupFctr = "d", breaks = 5, 
                   n.sims = 300, level = 0.9)
#> Warning: executing %dopar% sequentially: no parallel backend registered
impSim
#>   case bin   AvgFit     AvgFitSE nobs
#> 1    1   1 2.767921 2.847330e-04  193
#> 2    1   2 3.253321 6.623160e-05  240
#> 3    1   3 3.544864 5.988588e-05  254
#> 4    1   4 3.837757 7.343610e-05  265
#> 5    1   5 4.205953 1.836490e-04  176
```

The result of `REimpact` shows the change in the `yhat` as the case we
supplied to `newdata` is moved from the first to the fifth quintile in
terms of the magnitude of the group factor coefficient. We can see here
that the individual professor effect has a strong impact on the outcome
variable. This can be shown graphically as well:

``` r
ggplot(impSim, aes(x = factor(bin), y = AvgFit, ymin = AvgFit - 1.96*AvgFitSE, 
                   ymax = AvgFit + 1.96*AvgFitSE)) + 
  geom_pointrange() + theme_bw() + labs(x = "Bin of `d` term", y = "Predicted Fit")
```

![](man/figures/README_reImpactplot-1.png)<!-- -->

Here the standard error is a bit different – it is the weighted standard
error of the mean effect within the bin. It does not take into account
the variability within the effects of each observation in the bin –
accounting for this variation will be a future addition to `merTools`.

## Explore Substantive Impacts

Another feature of `merTools` is the ability to easily generate
hypothetical scenarios to explore the predicted outcomes of a `merMod`
object and understand what the model is saying in terms of the outcome
variable.

Let’s take the case where we want to explore the impact of a model with
an interaction term between a category and a continuous predictor.
First, we fit a model with interactions:

``` r
data(VerbAgg)
fmVA <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
           (1|id) + (1|item), family = binomial, 
           data = VerbAgg)
#> Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
#> Model failed to converge with max|grad| = 0.0543724 (tol = 0.002, component 1)
```

Now we prep the data using the `draw` function in `merTools`. Here we
draw the average observation from the model frame. We then `wiggle` the
data by expanding the dataframe to include the same observation repeated
but with different values of the variable specified by the `var`
parameter. Here, we expand the dataset to all values of `btype`, `situ`,
and `Anger` subsequently.

``` r
# Select the average case
newData <- draw(fmVA, type = "average")
newData <- wiggle(newData, varlist = "btype", 
                  valueslist = list(unique(VerbAgg$btype)))
newData <- wiggle(newData, var = "situ", 
                  valueslist = list(unique(VerbAgg$situ)))
newData <- wiggle(newData, var = "Anger", 
                  valueslist = list(unique(VerbAgg$Anger)))
head(newData, 10)
#>    r2 Anger Gender btype  situ id        item
#> 1   N    20      F curse other  5 S3WantCurse
#> 2   N    20      F scold other  5 S3WantCurse
#> 3   N    20      F shout other  5 S3WantCurse
#> 4   N    20      F curse  self  5 S3WantCurse
#> 5   N    20      F scold  self  5 S3WantCurse
#> 6   N    20      F shout  self  5 S3WantCurse
#> 7   N    11      F curse other  5 S3WantCurse
#> 8   N    11      F scold other  5 S3WantCurse
#> 9   N    11      F shout other  5 S3WantCurse
#> 10  N    11      F curse  self  5 S3WantCurse
```

The next step is familiar – we simply pass this new dataset to
`predictInterval` in order to generate predictions for these
counterfactuals. Then we plot the predicted values against the
continuous variable, `Anger`, and facet and group on the two categorical
variables `situ` and `btype` respectively.

``` r
plotdf <- predictInterval(fmVA, newdata = newData, type = "probability", 
            stat = "median", n.sims = 1000)
plotdf <- cbind(plotdf, newData)

ggplot(plotdf, aes(y = fit, x = Anger, color = btype, group = btype)) + 
  geom_point() + geom_smooth(aes(color = btype), method = "lm") + 
  facet_wrap(~situ) + theme_bw() +
  labs(y = "Predicted Probability")
#> `geom_smooth()` using formula 'y ~ x'
```

![](man/figures/README_substImpactPredict-1.png)<!-- -->

## Marginalizing Random Effects

``` r
# get cases
case_idx <- sample(1:nrow(VerbAgg), 10)
mfx <- REmargins(fmVA, newdata = VerbAgg[case_idx,], breaks = 4, groupFctr = "item", 
                 type = "probability")

ggplot(mfx, aes(y = fit_combined, x = breaks, group = case)) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  scale_y_continuous(breaks = 1:10/10, limits = c(0, 1)) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Quartile of item random effect Intercept for term 'item'", 
       y = "Predicted Probability", 
       title = "Simulated Effect of Item Intercept on Predicted Probability for 10 Random Cases")
```

![](man/figures/README_unnamed-chunk-15-1.png)<!-- -->
