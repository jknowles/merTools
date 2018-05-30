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
devtools::install_github("jknowles/merTools")
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
#>         fit      upr       lwr
#> 1  3.237114 5.178429 1.1223680
#> 2  3.044461 4.916916 1.0961706
#> 3  3.438357 5.332750 1.4219809
#> 4  3.027571 4.973256 0.9901124
#> 5  3.269985 5.462748 1.4136442
#> 6  3.286969 5.296241 1.3382832
#> 7  4.242488 6.542798 2.0041429
#> 8  3.721571 5.736674 1.8343516
#> 9  3.743598 5.891128 1.7195178
#> 10 3.419343 5.341034 1.4319983
```

Note that `predictInterval` is slower because it is computing simulations. It can also return all of the simulated `yhat` values as an attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily to draw the distributions of the parameters of the model. It then combines these simulated values to create a distribution of the `yhat` for each observation.

### Inspecting the Prediction Components

We can also explore the components of the prediction interval by asking `predictInterval` to return specific components of the prediction interval.

``` r
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 200, level = 0.9, 
                stat = 'median', which = "all")
#>      effect         fit      upr       lwr obs
#> 1  combined  3.13802196 4.719743  1.347545   1
#> 2  combined  2.95109536 4.872412  1.410548   2
#> 3  combined  3.37110677 5.469444  1.466620   3
#> 4  combined  3.16988494 5.195515  1.191728   4
#> 5  combined  3.40252170 5.463127  1.449450   5
#> 6  combined  3.15570458 5.303344  1.369328   6
#> 7  combined  4.22907656 6.020314  2.477142   7
#> 8  combined  3.77138257 5.908368  1.961115   8
#> 9  combined  3.67200743 5.702303  1.895175   9
#> 10 combined  3.39860415 5.780355  1.131460  10
#> 11        s  0.02263697 2.024355 -1.879137   1
#> 12        s  0.08910235 2.453711 -1.910101   2
#> 13        s  0.20390380 2.204549 -1.531280   3
#> 14        s  0.07476931 2.209324 -1.635799   4
#> 15        s -0.27105952 1.597833 -1.960480   5
#> 16        s -0.14858409 2.373469 -2.071048   6
#> 17        s  0.36679666 2.167334 -1.799921   7
#> 18        s  0.11402859 2.031803 -1.733151   8
#> 19        s  0.47395479 2.392591 -1.389445   9
#> 20        s  0.28443184 2.295856 -1.528997  10
#> 21        d -0.14468649 2.031511 -2.079169   1
#> 22        d -0.09577183 1.728302 -2.223127   2
#> 23        d  0.10723798 2.159129 -1.833482   3
#> 24        d -0.12653441 1.794607 -2.164757   4
#> 25        d  0.15969618 2.228805 -2.052448   5
#> 26        d  0.04889610 1.936657 -1.930860   6
#> 27        d  0.68879320 2.383263 -1.000621   7
#> 28        d  0.20416870 2.021266 -1.587005   8
#> 29        d  0.11384126 2.124815 -1.704743   9
#> 30        d -0.25093957 1.666946 -2.176818  10
#> 31    fixed  3.15731377 5.161839  1.284808   1
#> 32    fixed  3.25919968 5.116450  1.447072   2
#> 33    fixed  3.28844467 5.029445  1.231928   3
#> 34    fixed  3.15963283 5.190806  1.037147   4
#> 35    fixed  3.33292708 4.782191  1.367022   5
#> 36    fixed  3.29658917 5.107461  1.308795   6
#> 37    fixed  3.25299127 5.139824  1.146397   7
#> 38    fixed  3.22639471 4.739550  1.491264   8
#> 39    fixed  3.28930903 5.200827  1.322727   9
#> 40    fixed  3.31087713 5.158709  1.673070  10
```

This can lead to some useful plotting:

``` r
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
#> 1 (Intercept)  3.22450167  3.22518005 0.01658874
#> 2    service1 -0.07035669 -0.07063585 0.01398532
#> 3   lectage.L -0.18492073 -0.18684245 0.01609972
#> 4   lectage.Q  0.02392426  0.02514560 0.01201339
#> 5   lectage.C -0.02489267 -0.02481988 0.01432019
#> 6   lectage^4 -0.02083579 -0.01913810 0.01351510
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
#> 1         s       1 (Intercept)  0.12267659  0.16448740 0.2921818
#> 2         s       2 (Intercept) -0.02701206 -0.05781139 0.3079535
#> 3         s       3 (Intercept)  0.27639010  0.27821389 0.3101216
#> 4         s       4 (Intercept)  0.31600959  0.29023980 0.3422365
#> 5         s       5 (Intercept)  0.08113060  0.08345528 0.3495597
#> 6         s       6 (Intercept)  0.07523517  0.06432870 0.2299852
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
#> 2         d          1 _Intercept  0.3944916 0.08665150 835.3004    74
#> 3         d          6 _Intercept -0.4428947 0.03901988 239.5364    21
#> 4         d          7 _Intercept  0.6562683 0.03717200 997.3570    88
#> 5         d          8 _Intercept -0.6430679 0.02210017 138.3445    12
#> 6         d         12 _Intercept  0.1902942 0.04024063 702.3412    62
#> 7         d         13 _Intercept  0.2497467 0.03216255 750.0176    66
```

A nice features `expectedRank` is that you can return the expected rank for all factors simultaneously and use them:

``` r
ranks <- expectedRank(m1)
head(ranks)
#>   groupFctr groupLevel       term    estimate  std.error       ER pctER
#> 2         s          1 _Intercept  0.16732718 0.08165628 1931.568    65
#> 3         s          2 _Intercept -0.04409513 0.09234202 1368.161    46
#> 4         s          3 _Intercept  0.30382117 0.05204066 2309.940    78
#> 5         s          4 _Intercept  0.24756074 0.06641674 2151.827    72
#> 6         s          5 _Intercept  0.05232307 0.08174092 1627.693    55
#> 7         s          6 _Intercept  0.10191618 0.06648369 1772.548    60

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
#> 1    1   1 2.778793 3.275142e-04  193
#> 2    1   2 3.240974 7.034638e-05  240
#> 3    1   3 3.521627 5.823202e-05  254
#> 4    1   4 3.810390 5.598813e-05  265
#> 5    1   5 4.187144 1.975556e-04  176
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
#> $checkConv, : Model failed to converge with max|grad| = 0.0534956 (tol =
#> 0.001, component 1)
```

Now we prep the data using the `draw` function in `merTools`. Here we draw the average observation from the model frame. We then `wiggle` the data by expanding the dataframe to include the same observation repeated but with different values of the variable specified by the `var` parameter. Here, we expand the dataset to all values of `btype`, `situ`, and `Anger` subsequently.

``` r
# Select the average case
newData <- draw(fmVA, type = "average")
newData <- wiggle(newData, var = "btype", values = unique(VerbAgg$btype))
newData <- wiggle(newData, var = "situ", values = unique(VerbAgg$situ))
newData <- wiggle(newData, var = "Anger", values = unique(VerbAgg$Anger))
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
