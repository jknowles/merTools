[![Travis-CI Build
Status](https://travis-ci.org/jknowles/merTools.png?branch=master)](https://travis-ci.org/jknowles/merTools)
[![Coverage
Status](https://coveralls.io/repos/jknowles/merTools/badge.svg?branch=master)](https://coveralls.io/r/jknowles/merTools?branch=master)
[![Github
Issues](http://githubbadges.herokuapp.com/jknowles/merTools/issues.svg)](https://github.com/jknowles/merTools/issues)
[![Pending
Pull-Requests](http://githubbadges.herokuapp.com/jknowles/merTools/pulls.svg?style=flat)](https://github.com/jknowles/merTools/pulls)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/merTools)](https://cran.r-project.org/package=merTools)
[![Downloads](http://cranlogs.r-pkg.org/badges/merTools)](https://cran.r-project.org/package=merTools)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# merTools

A package for getting the most of our multilevel models in R

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

### merTools 0.4.1

#### New Features

  - Standard errors reported by `merModList` functions now apply the
    Rubin correction for multiple imputation

#### Bug fixes

  - Contribution by Alex Whitworth (@alexWhitworth) adding error
    checking to plotting functions

### merTools 0.4.0

#### New Features

  - Added vignette on using multilevel models with multiply imputed data
  - Added `fixef` and `ranef` generics for `merModList` objects
  - Added `fastdisp` generic for `merModList`
  - Added `summary` generic for `merModList`
  - Added `print` generic for `merModList`
  - Documented all generics for `merModList` including examples and a
    new imputation vignette
  - Added `modelInfo` generic for `merMod` objects that provides simple
    summary stats about a whole model

#### Bug Fixes

  - Fix bug that returned NaN for `std.error` of a multiply imputed
    `merModList` when calling `modelRandEffStats`
  - Fixed bug in `REimpact` where some column names in `newdata` would
    prevent the prediction intervals from being computed correctly.
    Users will now be warned.
  - Fixed bug in `wiggle` where documentation incorrectly stated the
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
#> 3.146336 3.165211 3.398499 3.114248 3.320686 3.252670 4.180896 3.845218 
#>        9       10 
#> 3.779336 3.331012
```

With `predictInterval` we obtain predictions that are more like the
standard objects produced by `lm` and
`glm`:

``` r
#predictInterval(m1, newdata = InstEval[1:10, ]) # all other parameters are optional
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 500, level = 0.9, 
                stat = 'median')
#>         fit      upr       lwr
#> 1  3.131864 5.129219 1.0945339
#> 2  3.178221 5.148597 1.2658735
#> 3  3.379966 5.516817 1.4895417
#> 4  3.103167 5.081666 0.9233663
#> 5  3.277556 5.399948 1.2272530
#> 6  3.149349 5.333070 1.3867458
#> 7  4.260277 6.219553 2.1694261
#> 8  3.875295 5.814772 1.9408867
#> 9  3.696745 5.793889 1.7367258
#> 10 3.231097 5.418401 1.3991943
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
#>      effect         fit      upr       lwr obs
#> 1  combined  3.27426294 5.078645  1.404931   1
#> 2  combined  3.31181634 5.294470  1.133825   2
#> 3  combined  3.42991423 5.726454  1.372737   3
#> 4  combined  3.14106388 5.086958  1.139264   4
#> 5  combined  3.21851807 5.005576  1.613092   5
#> 6  combined  3.13385479 5.357638  1.479084   6
#> 7  combined  4.20473915 6.038695  2.131664   7
#> 8  combined  3.89199549 5.824092  1.742773   8
#> 9  combined  3.73120768 5.937519  1.796547   9
#> 10 combined  3.38347478 5.139712  1.438058  10
#> 11        s  0.27050123 2.179032 -1.695454   1
#> 12        s  0.24091511 2.057844 -1.849768   2
#> 13        s  0.23153007 2.109371 -2.159917   3
#> 14        s  0.17190688 2.268925 -1.838580   4
#> 15        s  0.03743789 1.896608 -2.265675   5
#> 16        s -0.18521618 1.898214 -1.942181   6
#> 17        s  0.22645227 1.979926 -1.630725   7
#> 18        s  0.31321219 2.488689 -1.600209   8
#> 19        s  0.22507632 2.358639 -1.708491   9
#> 20        s  0.24543798 1.980311 -1.901763  10
#> 21        d -0.19334941 1.840901 -2.124230   1
#> 22        d -0.28982541 1.955074 -2.204380   2
#> 23        d -0.01613680 2.010074 -1.785555   3
#> 24        d -0.15524168 1.798061 -2.164504   4
#> 25        d  0.04642275 1.693999 -2.015113   5
#> 26        d -0.01555388 1.691761 -1.933333   6
#> 27        d  0.41534161 2.583055 -1.475386   7
#> 28        d  0.08841543 1.958569 -1.505878   8
#> 29        d  0.15627408 2.213819 -1.680117   9
#> 30        d -0.24753404 1.628555 -2.032027  10
#> 31    fixed  3.00038250 5.345275  1.047073   1
#> 32    fixed  3.31010106 5.258032  1.406672   2
#> 33    fixed  3.14307498 4.984565  1.302140   3
#> 34    fixed  3.18420609 5.233990  1.246140   4
#> 35    fixed  3.20555941 4.843914  1.144443   5
#> 36    fixed  3.28471448 5.179227  1.397785   6
#> 37    fixed  3.35935343 5.724092  1.375902   7
#> 38    fixed  3.27071043 5.168950  1.554638   8
#> 39    fixed  3.32543843 4.851246  1.432455   9
#> 40    fixed  3.30121263 5.451804  1.619428  10
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

plotdf <- bind_rows(plotdf, plotdfb, .id = "residVar")
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
#> 1 (Intercept)  3.22354670  3.22486405 0.01908013
#> 2    service1 -0.06887981 -0.07062088 0.01308958
#> 3   lectage.L -0.18513370 -0.18833678 0.01918627
#> 4   lectage.Q  0.02389151  0.02342562 0.01117116
#> 5   lectage.C -0.02465144 -0.02565338 0.01389634
#> 6   lectage^4 -0.02051981 -0.01880618 0.01178279
```

And we can also plot
this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](man/figures/README_FEsimPlot-1.png)<!-- -->

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.13373397  0.14351276 0.3554112
#> 2         s       2 (Intercept) -0.04046529 -0.01533224 0.3009859
#> 3         s       3 (Intercept)  0.25991981  0.25390020 0.3045212
#> 4         s       4 (Intercept)  0.21929110  0.15633635 0.3376836
#> 5         s       5 (Intercept)  0.07437984  0.06909753 0.3154968
#> 6         s       6 (Intercept)  0.09368613  0.12836103 0.2201394
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
#>   groupFctr groupLevel       term   estimate  std.error       ER pctER
#> 2         d          1 _Intercept  0.3944916 0.08665149 835.3004    74
#> 3         d          6 _Intercept -0.4428947 0.03901987 239.5364    21
#> 4         d          7 _Intercept  0.6562683 0.03717200 997.3570    88
#> 5         d          8 _Intercept -0.6430679 0.02210017 138.3445    12
#> 6         d         12 _Intercept  0.1902942 0.04024063 702.3412    62
#> 7         d         13 _Intercept  0.2497467 0.03216254 750.0176    66
```

A nice features `expectedRank` is that you can return the expected rank
for all factors simultaneously and use them:

``` r
ranks <- expectedRank(m1)
head(ranks)
#>   groupFctr groupLevel       term    estimate  std.error       ER pctER
#> 2         s          1 _Intercept  0.16732725 0.08165631 1931.569    65
#> 3         s          2 _Intercept -0.04409515 0.09234206 1368.161    46
#> 4         s          3 _Intercept  0.30382125 0.05204068 2309.940    78
#> 5         s          4 _Intercept  0.24756083 0.06641676 2151.827    72
#> 6         s          5 _Intercept  0.05232309 0.08174095 1627.693    55
#> 7         s          6 _Intercept  0.10191622 0.06648371 1772.548    60

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
#> 1    1   1 2.796121 2.839438e-04  193
#> 2    1   2 3.262527 7.421240e-05  240
#> 3    1   3 3.553761 5.812397e-05  254
#> 4    1   4 3.846521 6.754321e-05  265
#> 5    1   5 4.226088 1.763781e-04  176
```

The result of `REimpact` shows the change in the `yhat` as the case we
supplied to `newdata` is moved from the first to the fifth quintile in
terms of the magnitude of the group factor coefficient. We can see here
that the individual professor effect has a strong impact on the outcome
variable. This can be shown graphically as
well:

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
#> Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
#> $checkConv, : Model failed to converge with max|grad| = 0.0505464 (tol =
#> 0.001, component 1)
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
#>    r2 Anger Gender btype  situ  id        item
#> 1   N    20      F curse other 149 S3WantCurse
#> 2   N    20      F scold other 149 S3WantCurse
#> 3   N    20      F shout other 149 S3WantCurse
#> 4   N    20      F curse  self 149 S3WantCurse
#> 5   N    20      F scold  self 149 S3WantCurse
#> 6   N    20      F shout  self 149 S3WantCurse
#> 7   N    11      F curse other 149 S3WantCurse
#> 8   N    11      F scold other 149 S3WantCurse
#> 9   N    11      F shout other 149 S3WantCurse
#> 10  N    11      F curse  self 149 S3WantCurse
```

The next step is familiar – we simply pass this new dataset to
`predictInterval` in order to generate predictions for these
counterfactuals. Then we plot the predicted values against the
continuous variable, `Anger`, and facet and group on the two categorical
variables `situ` and `btype`
respectively.

``` r
plotdf <- predictInterval(fmVA, newdata = newData, type = "probability", 
            stat = "median", n.sims = 1000)
plotdf <- cbind(plotdf, newData)

ggplot(plotdf, aes(y = fit, x = Anger, color = btype, group = btype)) + 
  geom_point() + geom_smooth(aes(color = btype), method = "lm") + 
  facet_wrap(~situ) + theme_bw() +
  labs(y = "Predicted Probability")
```

![](man/figures/README_substImpactPredict-1.png)<!-- -->
