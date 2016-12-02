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

-   Rebuilt the readme.md to include new information about new features
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

### merTools 0.2.1

-   Fix dependency on `lme4` to ensure compatibility with latest changes.

### Version 0.2

#### New Functionality

-   Substantial performance enhancement for `predictInterval` which includes better handling of large numbers of parameters and simulations, performance tweaks for added speed (~10x), and parallel backend support (currently not optimized)
-   Add support for `probit` models and limited support for other `glmm` link functions, with warning (still do not know how to handle sigma parameter for these)
-   Add ability for user-specified seed for reproducibility
-   Add support for `blmerMod` and `bglmerMod` objects from the `blme` package
-   Add a `merModList` object for lists of `merMod` objects fitted to subsets of a dataset, useful for imputation or for working with extremely large datasets
-   Add a `print` method for `merModList` to mimic output of `summary.merMod`
-   Add a `VarCorr` method for `merModList`
-   Add new package data to demonstrate replication from selected published texts on multilevel modeling using different software (1982 High School and Beyond Survey data)

#### Other changes

-   Changed the default `n.sims` for the `predictInterval` function from 100 to 1,000 to give better coverage and reflect performance increase
-   Changed the default for `level` in `predictInterval` to be 0.8 instead of 0.95 to reflect that 0.95 prediction intervals are more conservative than most users need

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

![](readmeplot/README-predPanel.png)

On the first tab, the function presents the prediction intervals for the data selected by user which are calculated using the `predictInterval` function within the package. This function calculates prediction intervals quickly by sampling from the simulated distribution of the fixed effect and random effect terms and combining these simulated estimates to produce a distribution of predictions for each observation. This allows prediction intervals to be generated from very large models where the use of `bootMer` would not be feasible computationally.

![](readmeplot/README-effPanel.png)

On the next tab the distribution of the fixed effect and group-level effects is depicted on confidence interval plots. These are useful for diagnostics and provide a way to inspect the relative magnitudes of various parameters. This tab makes use of four related functions in `merTools`: `FEsim`, `plotFEsim`, `REsim` and `plotREsim` which are available to be used on their own as well.

![](readmeplot/README-substPanel.png)

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
#> 1  3.132867 5.157970 1.183171
#> 2  3.180017 5.072290 0.973230
#> 3  3.439371 5.325379 1.358078
#> 4  3.131259 5.142101 1.236331
#> 5  3.426534 5.172656 1.325485
#> 6  3.204123 5.279035 1.245494
#> 7  4.206803 6.138637 2.241488
#> 8  3.886083 5.982433 1.830014
#> 9  3.681801 5.725326 1.802468
#> 10 3.353775 5.504170 1.467024
```

Note that `predictInterval` is slower because it is computing simulations. It can also return all of the simulated `yhat` values as an attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily to draw the distributions of the parameters of the model. It then combines these simulated values to create a distribution of the `yhat` for each observation.

### Inspecting the Prediction Components

We can also explore the components of the prediction interval by asking `predictInterval` to return specific components of the prediction interval.

``` r
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 200, level = 0.9, 
                stat = 'median', which = "all")
#>      effect          fit      upr       lwr obs
#> 1  combined  3.242855583 5.256198  1.166928   1
#> 2  combined  3.281278745 5.327124  1.047201   2
#> 3  combined  3.540671486 5.477600  1.193292   3
#> 4  combined  3.115095384 5.172895  1.243576   4
#> 5  combined  3.108293629 5.241617  1.425906   5
#> 6  combined  3.341020492 5.019363  1.460654   6
#> 7  combined  4.225094994 5.749342  2.142919   7
#> 8  combined  3.768235942 6.134869  1.904707   8
#> 9  combined  3.780869784 5.688571  1.724311   9
#> 10 combined  3.270161364 5.617584  1.612083  10
#> 11        s  0.339715613 2.370676 -1.755473   1
#> 12        s  0.008750713 2.069383 -1.812534   2
#> 13        s  0.431764762 2.167529 -1.699329   3
#> 14        s  0.295036886 2.434658 -1.756879   4
#> 15        s -0.132865616 1.855530 -2.066101   5
#> 16        s -0.191921491 1.823088 -2.017646   6
#> 17        s  0.419851243 2.384942 -1.455827   7
#> 18        s  0.420740816 2.182718 -1.405249   8
#> 19        s  0.183052132 2.204805 -1.734104   9
#> 20        s  0.378577708 2.215289 -1.422765  10
#> 21        d -0.371815298 1.668305 -2.222380   1
#> 22        d -0.256243810 1.791162 -2.021152   2
#> 23        d  0.009742000 2.166273 -1.891261   3
#> 24        d -0.363916557 1.463719 -2.132313   4
#> 25        d  0.076165324 2.162862 -1.973062   5
#> 26        d -0.057323390 2.092815 -1.694397   6
#> 27        d  0.787872063 2.697186 -1.218455   7
#> 28        d  0.303469971 2.269464 -1.591569   8
#> 29        d  0.225490500 1.793513 -1.902528   9
#> 30        d -0.361997595 1.462713 -2.385442  10
#> 31    fixed  3.197006848 5.183858  1.231755   1
#> 32    fixed  3.182668601 4.967376  1.691922   2
#> 33    fixed  3.337971033 5.462194  1.328967   3
#> 34    fixed  3.003346414 5.144884  1.075642   4
#> 35    fixed  3.284434481 5.086225  1.539723   5
#> 36    fixed  3.406863154 4.972823  1.415661   6
#> 37    fixed  3.161917368 5.388678  1.031420   7
#> 38    fixed  3.303104960 5.137719  1.314559   8
#> 39    fixed  3.259537665 4.996122  1.426116   9
#> 40    fixed  3.173104292 5.081483  1.196755  10
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

![](readmeplot/README-unnamed-chunk-8-1.png)

We can also investigate the makeup of the prediction for each observation.

``` r
ggplot(plotdf[plotdf$obs < 6,], 
       aes(x = effect, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  facet_grid(residVar~obs) + theme_bw()
```

![](readmeplot/README-unnamed-chunk-9-1.png)

Plotting
--------

`merTools` also provides functionality for inspecting `merMod` objects visually. The easiest are getting the posterior distributions of both fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22172541  3.22127980 0.01747097
#> 2    service1 -0.07056833 -0.07136218 0.01531774
#> 3   lectage.L -0.18996987 -0.19067454 0.01428177
#> 4   lectage.Q  0.02393843  0.02261053 0.01173272
#> 5   lectage.C -0.02527826 -0.02578302 0.01416652
#> 6   lectage^4 -0.02031094 -0.02174745 0.01122081
```

And we can also plot this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](readmeplot/README-FEsimPlot-1.png)

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.22599246  0.25097098 0.2985394
#> 2         s       2 (Intercept) -0.06229065 -0.11287415 0.2851206
#> 3         s       3 (Intercept)  0.26950589  0.25761638 0.2817762
#> 4         s       4 (Intercept)  0.22514387  0.21662447 0.2866336
#> 5         s       5 (Intercept)  0.06392735  0.08611845 0.3112211
#> 6         s       6 (Intercept)  0.09866372  0.08245427 0.2644501
```

``` r
plotREsim(REsim(m1, n.sims = 100), stat = 'median', sd = TRUE)
```

![](readmeplot/README-reSimplot-1.png)

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
#> 2         s          1 _Intercept  0.16732717 0.08165627 1931.568    65
#> 3         s          2 _Intercept -0.04409512 0.09234201 1368.161    46
#> 4         s          3 _Intercept  0.30382116 0.05204066 2309.940    78
#> 5         s          4 _Intercept  0.24756072 0.06641674 2151.827    72
#> 6         s          5 _Intercept  0.05232306 0.08174091 1627.693    55
#> 7         s          6 _Intercept  0.10191617 0.06648369 1772.548    60

ggplot(ranks, aes(x = term, y = estimate)) + 
  geom_violin(fill = "gray50") + facet_wrap(~groupFctr) +
  theme_bw()
```

![](readmeplot/README-unnamed-chunk-13-1.png)

Effect Simulation
-----------------

It can still be difficult to interpret the results of LMM and GLMM models, especially the relative influence of varying parameters on the predicted outcome. This is where the `REimpact` and the `wiggle` functions in `merTools` can be handy.

``` r
impSim <- REimpact(m1, InstEval[7, ], groupFctr = "d", breaks = 5, 
                   n.sims = 300, level = 0.9)
#> Warning: executing %dopar% sequentially: no parallel backend registered
impSim
#>   case bin   AvgFit     AvgFitSE nobs
#> 1    1   1 2.791733 3.178347e-04  193
#> 2    1   2 3.268277 7.266088e-05  240
#> 3    1   3 3.565191 5.879801e-05  254
#> 4    1   4 3.850698 6.746706e-05  265
#> 5    1   5 4.238028 1.837003e-04  176
```

The result of `REimpact` shows the change in the `yhat` as the case we supplied to `newdata` is moved from the first to the fifth quintile in terms of the magnitude of the group factor coefficient. We can see here that the individual professor effect has a strong impact on the outcome variable. This can be shown graphically as well:

``` r
ggplot(impSim, aes(x = factor(bin), y = AvgFit, ymin = AvgFit - 1.96*AvgFitSE, 
                   ymax = AvgFit + 1.96*AvgFitSE)) + 
  geom_pointrange() + theme_bw() + labs(x = "Bin of `d` term", y = "Predicted Fit")
```

![](readmeplot/README-reImpactplot-1.png)

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
#> $checkConv, : Model failed to converge with max|grad| = 0.028633 (tol =
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
#> 1   N    20      F curse other  3 S3WantCurse
#> 2   N    20      F scold other  3 S3WantCurse
#> 3   N    20      F shout other  3 S3WantCurse
#> 4   N    20      F curse  self  3 S3WantCurse
#> 5   N    20      F scold  self  3 S3WantCurse
#> 6   N    20      F shout  self  3 S3WantCurse
#> 7   N    11      F curse other  3 S3WantCurse
#> 8   N    11      F scold other  3 S3WantCurse
#> 9   N    11      F shout other  3 S3WantCurse
#> 10  N    11      F curse  self  3 S3WantCurse
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

![](readmeplot/README-substImpactPredict-1.png)
