<!-- README.md is generated from README.Rmd. Please edit that file -->
merTools
========

A package for getting the most of our multilevel models in R

by Jared E. Knowles and Carl Frederick

[![Travis-CI Build Status](https://travis-ci.org/jknowles/merTools.png?branch=master)](https://travis-ci.org/jknowles/merTools) [![Coverage Status](https://coveralls.io/repos/jknowles/merTools/badge.svg?branch=master)](https://coveralls.io/r/jknowles/merTools?branch=master) [![Github Issues](http://githubbadges.herokuapp.com/jknowles/merTools/issues.svg)](https://github.com/jknowles/merTools/issues) [![Pending Pull-Requests](http://githubbadges.herokuapp.com/jknowles/merTools/pulls.svg?style=flat)](https://github.com/jknowles/merTools/pulls) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/merTools)](https://cran.r-project.org/package=merTools) [![Downloads](http://cranlogs.r-pkg.org/badges/merTools)](https://cran.r-project.org/package=merTools)

Working with generalized linear mixed models (GLMM) and linear mixed models (LMM) has become increasingly easy with advances in the `lme4` package. As we have found ourselves using these models more and more within our work, we, the authors, have developed a set of tools for simplifying and speeding up common tasks for interacting with `merMod` objects from `lme4`. This package provides those tools.

Installation
------------

``` r
# development version
library(devtools)
install_github("jknowles/merTools")

# CRAN version
install.packages("merTools")
```

Recent Updates
--------------

### merTools 0.3.0

-   Improve handling of formulas. If the original `merMod` has functions specified in the formula, the `draw` and `wiggle` functions will check for this and attempt to respect these variable transformations. Where this is not possible a warning will be issued. Most common transformations are respected as long as the the original variable is passed untransformed to the model.
-   Change the calculations of the residual variance. Previously residual variance was used to inflate both the variance around the fixed parameters and around the predicted values themselves. This was incorrect and resulted in overly conservative estimates. Now the residual variance is appropriately only used around the final predictions
-   New option for `predictInterval` that allows the user to return the full interval, the fixed component, the random component, or the fixed and each random component separately for each observation
-   Fixed a bug with slope+intercept random terms that caused a miscalculation of the random component
-   Add comparison to `rstanarm` to the Vignette
-   Make `expectedRank` output more `tidy` like and allow function to calculate expected rank for all terms at once
-   Note, this breaks the API by changing the names of the columns in the output of this function
-   Remove tests that test for timing to avoid issues with R-devel JIT compiler
-   Remove `plyr` and replace with `dplyr`
-   Fix issue \#62 `varList` will now throw an error if `==` is used instead of `=`
-   Fix issue \#54 `predictInterval` did not included random effects in calculations when `newdata` had more than 1000 rows and/or user specified `parallel=TRUE`. Note: fix was to disable the `.paropts` option for `predictInterval` ... user can still specify for *temporary* backward compatibility but this should be either removed or fixed in the permanent solution.
-   Fix issue \#53 about problems with `predictInterval` when only specific levels of a grouping factor are in `newdata` with the colon specification of interactions
-   Fix issue \#52 ICC wrong calculations ... we just needed to square the standard deviations that we pulled

See [NEWS.md](https://github.com/jknowles/merTools/blob/master/NEWS.md) for more details.

Shiny App and Demo
------------------

The easiest way to demo the features of this application is to use the bundled Shiny application which launches a number of the metrics here to aide in exploring the model. To do this:

``` r
library(merTools)
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data
```

![](tools/readmeplot/README-predPanel.png)

On the first tab, the function presents the prediction intervals for the data selected by user which are calculated using the `predictInterval` function within the package. This function calculates prediction intervals quickly by sampling from the simulated distribution of the fixed effect and random effect terms and combining these simulated estimates to produce a distribution of predictions for each observation. This allows prediction intervals to be generated from very large models where the use of `bootMer` would not be feasible computationally.

![](tools/readmeplot/README-effPanel.png)

On the next tab the distribution of the fixed effect and group-level effects is depicted on confidence interval plots. These are useful for diagnostics and provide a way to inspect the relative magnitudes of various parameters. This tab makes use of four related functions in `merTools`: `FEsim`, `plotFEsim`, `REsim` and `plotREsim` which are available to be used on their own as well.

![](tools/readmeplot/README-substPanel.png)

On the third tab are some convenient ways to show the influence or magnitude of effects by leveraging the power of `predictInterval`. For each case, up to 12, in the selected data type, the user can view the impact of changing either one of the fixed effect or one of the grouping level terms. Using the `REimpact` function, each case is simulated with the model's prediction if all else was held equal, but the observation was moved through the distribution of the fixed effect or the random effect term. This is plotted on the scale of the dependent variable, which allows the user to compare the magnitude of effects across variables, and also between models on the same data.

Predicting
----------

Standard prediction looks like so.

``` r
predict(m1, newdata = InstEval[1:10, ])
#>        1        2        3        4        5        6        7        8 
#> 3.146336 3.165211 3.398499 3.114248 3.320686 3.252670 4.180896 3.845218 
#>        9       10 
#> 3.779336 3.331012
```

With `predictInterval` we obtain predictions that are more like the standard objects produced by `lm` and `glm`:

``` r
#predictInterval(m1, newdata = InstEval[1:10, ]) # all other parameters are optional
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 500, level = 0.9, 
                stat = 'median')
#>         fit      upr      lwr
#> 1  3.168065 5.101014 1.158220
#> 2  3.145478 5.274931 1.190373
#> 3  3.392613 5.590413 1.539565
#> 4  3.138387 5.028683 1.311576
#> 5  3.397307 5.486894 1.287190
#> 6  3.263650 5.191198 1.199061
#> 7  4.251013 6.065739 2.107230
#> 8  4.047291 5.813104 1.974133
#> 9  3.686001 5.838891 1.816260
#> 10 3.264626 5.234812 1.289614
```

Note that `predictInterval` is slower because it is computing simulations. It can also return all of the simulated `yhat` values as an attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily to draw the distributions of the parameters of the model. It then combines these simulated values to create a distribution of the `yhat` for each observation.

### Inspecting the Prediction Components

We can also explore the components of the prediction interval by asking `predictInterval` to return specific components of the prediction interval.

``` r
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 200, level = 0.9, 
                stat = 'median', which = "all")
#>      effect          fit      upr       lwr obs
#> 1  combined  3.112867407 4.986883  1.041417   1
#> 2  combined  2.989389447 4.991738  1.142979   2
#> 3  combined  3.418377391 5.446396  1.415533   3
#> 4  combined  3.204571267 5.170244  1.141322   4
#> 5  combined  3.285116884 5.267212  1.540779   5
#> 6  combined  3.082323400 5.250609  0.898447   6
#> 7  combined  4.169363868 6.190107  2.240950   7
#> 8  combined  3.827953300 6.013183  1.708390   8
#> 9  combined  3.708192514 5.808735  1.784811   9
#> 10 combined  3.379363265 5.358377  1.147798  10
#> 11        s  0.214052526 1.878793 -1.943308   1
#> 12        s  0.200860375 1.844383 -1.649367   2
#> 13        s  0.130458391 2.324902 -1.656628   3
#> 14        s  0.112916598 2.196385 -1.855981   4
#> 15        s -0.147686282 1.784524 -1.789859   5
#> 16        s  0.079661473 1.831603 -1.902383   6
#> 17        s  0.348544074 2.183217 -1.909237   7
#> 18        s  0.415829437 2.239591 -1.603185   8
#> 19        s  0.205510317 2.134001 -1.900833   9
#> 20        s  0.227854238 2.225521 -1.532491  10
#> 21        d -0.345513817 1.954967 -2.021580   1
#> 22        d -0.209006946 1.924467 -2.172302   2
#> 23        d  0.177014731 1.922807 -1.728790   3
#> 24        d -0.092382741 1.743827 -2.195150   4
#> 25        d  0.093043286 1.902361 -1.887277   5
#> 26        d  0.009651465 2.064675 -1.840032   6
#> 27        d  0.756918178 2.336602 -1.092124   7
#> 28        d  0.344379468 2.277583 -1.497924   8
#> 29        d  0.226951343 1.991785 -1.799554   9
#> 30        d -0.206121277 1.887244 -2.098502  10
#> 31    fixed  3.285346958 5.349092  1.512967   1
#> 32    fixed  3.414238014 5.085391  1.402841   2
#> 33    fixed  3.191911826 5.025254  1.166261   3
#> 34    fixed  3.054363384 4.911773  1.128827   4
#> 35    fixed  3.306328559 5.134675  1.692722   5
#> 36    fixed  3.258911806 5.152391  1.678185   6
#> 37    fixed  2.974008141 4.851762  1.068538   7
#> 38    fixed  3.359217613 5.265864  1.473522   8
#> 39    fixed  3.124721362 4.933596  1.376023   9
#> 40    fixed  3.309916523 5.050473  1.252730  10
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

![](tools/readmeplot/README-unnamed-chunk-8-1.png)

We can also investigate the makeup of the prediction for each observation.

``` r
ggplot(plotdf[plotdf$obs < 6,], 
       aes(x = effect, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  facet_grid(residVar~obs) + theme_bw()
```

![](tools/readmeplot/README-unnamed-chunk-9-1.png)

Plotting
--------

`merTools` also provides functionality for inspecting `merMod` objects visually. The easiest are getting the posterior distributions of both fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22774912  3.22942908 0.02039429
#> 2    service1 -0.07170745 -0.07213702 0.01163877
#> 3   lectage.L -0.18569377 -0.18279299 0.01565488
#> 4   lectage.Q  0.02178664  0.02132772 0.01164987
#> 5   lectage.C -0.02577289 -0.02712761 0.01257388
#> 6   lectage^4 -0.02243780 -0.02179817 0.01154704
```

And we can also plot this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](tools/readmeplot/README-FEsimPlot-1.png)

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.10840825  0.12737177 0.3052798
#> 2         s       2 (Intercept) -0.08274821 -0.09812382 0.3369063
#> 3         s       3 (Intercept)  0.29195086  0.33069702 0.2673189
#> 4         s       4 (Intercept)  0.23319286  0.26716723 0.3119964
#> 5         s       5 (Intercept)  0.05002205  0.02259669 0.3061090
#> 6         s       6 (Intercept)  0.12910788  0.11464799 0.2531423
```

``` r
plotREsim(REsim(m1, n.sims = 100), stat = 'median', sd = TRUE)
```

![](tools/readmeplot/README-reSimplot-1.png)

Note that `plotREsim` highlights group levels that have a simulated distribution that does not overlap 0 -- these appear darker. The lighter bars represent grouping levels that are not distinguishable from 0 in the data.

Sometimes the random effects can be hard to interpret and not all of them are meaningfully different from zero. To help with this `merTools` provides the `expectedRank` function, which provides the percentile ranks for the observed groups in the random effect distribution taking into account both the magnitude and uncertainty of the estimated effect for each group.

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

A nice features `expectedRank` is that you can return the expected rank for all factors simultaneously and use them:

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

![](tools/readmeplot/README-unnamed-chunk-13-1.png)

Effect Simulation
-----------------

It can still be difficult to interpret the results of LMM and GLMM models, especially the relative influence of varying parameters on the predicted outcome. This is where the `REimpact` and the `wiggle` functions in `merTools` can be handy.

``` r
impSim <- REimpact(m1, InstEval[7, ], groupFctr = "d", breaks = 5, 
                   n.sims = 300, level = 0.9)
#> Warning: executing %dopar% sequentially: no parallel backend registered
impSim
#>   case bin   AvgFit     AvgFitSE nobs
#> 1    1   1 2.784442 3.217872e-04  193
#> 2    1   2 3.257485 6.797772e-05  240
#> 3    1   3 3.556476 5.686517e-05  254
#> 4    1   4 3.838122 6.070734e-05  265
#> 5    1   5 4.217048 1.810544e-04  176
```

The result of `REimpact` shows the change in the `yhat` as the case we supplied to `newdata` is moved from the first to the fifth quintile in terms of the magnitude of the group factor coefficient. We can see here that the individual professor effect has a strong impact on the outcome variable. This can be shown graphically as well:

``` r
ggplot(impSim, aes(x = factor(bin), y = AvgFit, ymin = AvgFit - 1.96*AvgFitSE, 
                   ymax = AvgFit + 1.96*AvgFitSE)) + 
  geom_pointrange() + theme_bw() + labs(x = "Bin of `d` term", y = "Predicted Fit")
```

![](tools/readmeplot/README-reImpactplot-1.png)

Here the standard error is a bit different -- it is the weighted standard error of the mean effect within the bin. It does not take into account the variability within the effects of each observation in the bin -- accounting for this variation will be a future addition to `merTools`.

Explore Substantive Impacts
---------------------------

Another feature of `merTools` is the ability to easily generate hypothetical scenarios to explore the predicted outcomes of a `merMod` object and understand what the model is saying in terms of the outcome variable.

Let's take the case where we want to explore the impact of a model with an interaction term between a category and a continuous predictor. First, we fit a model with interactions:

``` r
data(VerbAgg)
fmVA <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
           (1|id) + (1|item), family = binomial, 
           data = VerbAgg)
#> Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
#> $checkConv, : Model failed to converge with max|grad| = 0.0505464 (tol =
#> 0.001, component 1)
```

Now we prep the data using the `draw` function in `merTools`. Here we draw the average observation from the model frame. We then `wiggle` the data by expanding the dataframe to include the same observation repeated but with different values of the variable specified by the `var` parameter. Here, we expand the dataset to all values of `btype`, `situ`, and `Anger` subsequently.

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

The next step is familiar -- we simply pass this new dataset to `predictInterval` in order to generate predictions for these counterfactuals. Then we plot the predicted values against the continuous variable, `Anger`, and facet and group on the two categorical variables `situ` and `btype` respectively.

``` r
plotdf <- predictInterval(fmVA, newdata = newData, type = "probability", 
            stat = "median", n.sims = 1000)
plotdf <- cbind(plotdf, newData)

ggplot(plotdf, aes(y = fit, x = Anger, color = btype, group = btype)) + 
  geom_point() + geom_smooth(aes(color = btype), method = "lm") + 
  facet_wrap(~situ) + theme_bw() +
  labs(y = "Predicted Probability")
```

![](tools/readmeplot/README-substImpactPredict-1.png)
