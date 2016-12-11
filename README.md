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
#> 1  3.114604 5.156392 1.081717
#> 2  3.117272 5.030120 1.242363
#> 3  3.423697 5.532564 1.322712
#> 4  3.083484 5.122665 1.188285
#> 5  3.213532 5.172851 1.192604
#> 6  3.096635 5.165586 1.265317
#> 7  4.154033 6.207856 2.062457
#> 8  3.874819 5.947561 1.821840
#> 9  3.819173 5.881570 1.796024
#> 10 3.255822 5.234510 1.348180
```

Note that `predictInterval` is slower because it is computing simulations. It can also return all of the simulated `yhat` values as an attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily to draw the distributions of the parameters of the model. It then combines these simulated values to create a distribution of the `yhat` for each observation.

### Inspecting the Prediction Components

We can also explore the components of the prediction interval by asking `predictInterval` to return specific components of the prediction interval.

``` r
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 200, level = 0.9, 
                stat = 'median', which = "all")
#>      effect         fit      upr        lwr obs
#> 1  combined  2.94585696 5.021224  0.9469043   1
#> 2  combined  3.27155895 4.976318  1.4859477   2
#> 3  combined  3.47215103 5.482550  1.4411054   3
#> 4  combined  3.20818135 5.115176  1.2518959   4
#> 5  combined  3.30141732 5.125936  1.4969704   5
#> 6  combined  3.18741956 5.012293  1.1022775   6
#> 7  combined  4.02052064 6.112281  2.1401665   7
#> 8  combined  3.81252578 5.970377  1.6762839   8
#> 9  combined  3.95422686 5.724368  1.8217968   9
#> 10 combined  3.37779081 5.589283  1.4468875  10
#> 11        s  0.15184349 2.192223 -1.5854914   1
#> 12        s  0.18125006 2.078097 -1.8952277   2
#> 13        s  0.11049947 2.857706 -1.9754246   3
#> 14        s  0.18872857 2.182267 -1.9760269   4
#> 15        s  0.01356864 1.855524 -2.0058637   5
#> 16        s -0.08246566 1.552007 -1.9481578   6
#> 17        s  0.45712964 2.365747 -1.6018123   7
#> 18        s  0.35318191 2.466022 -1.7020094   8
#> 19        s -0.04810836 2.443321 -1.6521329   9
#> 20        s  0.48134347 2.294846 -1.5030829  10
#> 21        d -0.23909115 1.953999 -2.1678440   1
#> 22        d -0.18816371 1.894217 -2.1977454   2
#> 23        d -0.18722419 1.947029 -2.0175046   3
#> 24        d -0.29574341 1.689394 -2.0905452   4
#> 25        d  0.00607903 1.955012 -1.8320733   5
#> 26        d  0.08034936 1.767552 -1.4940104   6
#> 27        d  0.72769488 2.719575 -1.1363420   7
#> 28        d  0.19514852 1.934265 -1.6582669   8
#> 29        d  0.22215313 2.144313 -1.9546232   9
#> 30        d -0.32046741 1.241536 -2.4551447  10
#> 31    fixed  3.15255529 5.103040  1.4456085   1
#> 32    fixed  3.15208224 4.971556  1.5238647   2
#> 33    fixed  3.07414134 4.965873  1.3234276   3
#> 34    fixed  3.15997183 4.854703  1.1075294   4
#> 35    fixed  3.31791081 5.248444  1.3468899   5
#> 36    fixed  3.28681289 5.262730  1.0778400   6
#> 37    fixed  3.21353040 5.345447  1.4695929   7
#> 38    fixed  3.37980340 5.339094  1.4238363   8
#> 39    fixed  3.34750489 5.124405  1.4098634   9
#> 40    fixed  3.14722350 5.331178  1.2528532  10
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
#> 1 (Intercept)  3.22451072  3.22467380 0.01986621
#> 2    service1 -0.07276963 -0.07311372 0.01325825
#> 3   lectage.L -0.18689581 -0.18670767 0.01655902
#> 4   lectage.Q  0.02340434  0.02353874 0.01168251
#> 5   lectage.C -0.02371317 -0.02476775 0.01324309
#> 6   lectage^4 -0.01940582 -0.01863110 0.01268280
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
#> 1         s       1 (Intercept)  0.15588013  0.20416978 0.3020464
#> 2         s       2 (Intercept) -0.05444166 -0.05035208 0.3477367
#> 3         s       3 (Intercept)  0.33650172  0.35410975 0.2983894
#> 4         s       4 (Intercept)  0.24081609  0.22966938 0.2836278
#> 5         s       5 (Intercept)  0.04425092  0.04415465 0.3224734
#> 6         s       6 (Intercept)  0.09972308  0.08947532 0.2295493
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
#> 2         d          1 _Intercept  0.3944915 0.08665148 835.3004    74
#> 3         d          6 _Intercept -0.4428947 0.03901987 239.5364    21
#> 4         d          7 _Intercept  0.6562682 0.03717199 997.3570    88
#> 5         d          8 _Intercept -0.6430679 0.02210017 138.3444    12
#> 6         d         12 _Intercept  0.1902942 0.04024062 702.3412    62
#> 7         d         13 _Intercept  0.2497466 0.03216254 750.0176    66
```

A nice features `expectedRank` is that you can return the expected rank for all factors simultaneously and use them:

``` r
ranks <- expectedRank(m1)
head(ranks)
#>   groupFctr groupLevel       term    estimate  std.error       ER pctER
#> 2         s          1 _Intercept  0.16732726 0.08165631 1931.569    65
#> 3         s          2 _Intercept -0.04409515 0.09234207 1368.161    46
#> 4         s          3 _Intercept  0.30382126 0.05204068 2309.940    78
#> 5         s          4 _Intercept  0.24756085 0.06641676 2151.827    72
#> 6         s          5 _Intercept  0.05232309 0.08174096 1627.693    55
#> 7         s          6 _Intercept  0.10191623 0.06648372 1772.548    60

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
#> 1    1   1 2.771751 2.882702e-04  193
#> 2    1   2 3.241434 7.178053e-05  240
#> 3    1   3 3.529431 5.378721e-05  254
#> 4    1   4 3.816805 6.711257e-05  265
#> 5    1   5 4.187946 2.075199e-04  176
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
