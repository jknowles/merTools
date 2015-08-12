<!-- README.md is generated from README.Rmd. Please edit that file -->
merTools
========

A package for getting the most of our multilevel models in R

by Jared E. Knowles and Carl Frederick

[![Travis-CI Build Status](https://travis-ci.org/jknowles/merTools.png?branch=master)](https://travis-ci.org/jknowles/merTools) [![Coverage Status](https://coveralls.io/repos/jknowles/merTools/badge.svg?branch=master)](https://coveralls.io/r/jknowles/merTools?branch=master) [![Github Issues](http://githubbadges.herokuapp.com/jknowles/merTools/issues.svg)](https://github.com/jknowles/merTools/issues) [![Pending Pull-Requests](http://githubbadges.herokuapp.com/jknowles/merTools/pulls.svg?style=flat)](https://github.com/jknowles/merTools/pulls) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/caretEnsemble)](http://cran.r-project.org/web/packages/merTools) [![Downloads](http://cranlogs.r-pkg.org/badges/merTools)](http://cran.rstudio.com/package=merTools)

Working with generalized linear mixed models (GLMM) and linear mixed models (LMM) has become increasingly easy with advances in the `lme4` package. As we have found ourselves using these models more and more within our work, we, the authors, have developed a set of tools for simplifying and speeding up common tasks for interacting with `merMod` objects from `lme4`. This package provides those tools.

Installation
------------

``` r
# development version
library(devtools)
install_github("jknowles/merTools")

# CRAN version -- coming soon
install.packages("merTools")
```

Shiny App and Demo
------------------

The easiest way to demo the features of this application is to use the bundled Shiny application which launches a number of the metrics here to aide in exploring the model. To do this:

``` r
devtools::install_github("jknowles/merTools")
library(merTools)
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data
```

![](README-predPanel.png)

On the first tab, the function presents the prediction intervals for the data selected by user which are calculated using the `predictInterval` function within the package. This function calculates prediction intervals quickly by sampling from the simulated distribution of the fixed effect and random effect terms and combining these simulated estimates to produce a distribution of predictions for each observation. This allows prediction intervals to be generated from very large models where the use of `bootMer` would not be feasible computationally.

![](README-effPanel.png)

On the next tab the distribution of the fixed effect and group-level effects is depicted on confidence interval plots. These are useful for diagnostics and provide a way to inspect the relative magnitudes of various parameters. This tab makes use of four related functions in `merTools`: `FEsim`, `plotFEsim`, `REsim` and `plotREsim` which are available to be used on their own as well.

![](README-substPanel.png)

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
#>         fit      lwr      upr
#> 1  3.212678 1.063293 5.021020
#> 2  3.163909 1.438549 5.043916
#> 3  3.423594 1.450930 5.329387
#> 4  3.184908 1.082163 4.985627
#> 5  3.272646 1.538196 5.452890
#> 6  3.229198 1.295352 5.283061
#> 7  4.128929 2.228645 6.094850
#> 8  3.910574 1.648967 5.849607
#> 9  3.729399 1.889009 5.711374
#> 10 3.438017 1.418762 5.508910
```

Note that `predictInterval` is slower because it is computing simulations. It can also return all of the simulated `yhat` values as an attribute to the predict object itself.

`predictInterval` uses the `sim` function from the `arm` package heavily to draw the distributions of the parameters of the model. It then combines these simulated values to create a distribution of the `yhat` for each observation.

Plotting
--------

`merTools` also provides functionality for inspecting `merMod` objects visually. The easiest are getting the posterior distributions of both fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22102856  3.22094538 0.01857148
#> 2    service1 -0.06968367 -0.06987645 0.01362059
#> 3   lectage.L -0.18605240 -0.18445724 0.01570463
#> 4   lectage.Q  0.02434054  0.02300220 0.01334311
#> 5   lectage.C -0.02490739 -0.02606262 0.01192429
#> 6   lectage^4 -0.01891636 -0.01875975 0.01271187
```

And we can also plot this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](README-FEsimPlot-1.png)

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.22058081  0.21648847 0.3233950
#> 2         s       2 (Intercept) -0.10235168 -0.10254166 0.3149336
#> 3         s       3 (Intercept)  0.29383942  0.26797892 0.3291858
#> 4         s       4 (Intercept)  0.23623372  0.25423707 0.3007441
#> 5         s       5 (Intercept)  0.08476158  0.07277292 0.3131551
#> 6         s       6 (Intercept)  0.09851428  0.07116080 0.2581527
```

``` r
plotREsim(REsim(m1, n.sims = 100), stat = 'median', sd = TRUE)
```

![](README-reSimplot-1.png)

Note that `plotREsim` highlights group levels that have a simulated distribution that does not overlap 0 -- these appear darker. The lighter bars represent grouping levels that are not distinguishable from 0 in the data.

Sometimes the random effects can be hard to interpret and not all of them are meaningfully different from zero. To help with this `merTools` provides the `expectedRank` function, which provides the percentile ranks for the observed groups in the random effect distribution taking into account both the magnitude and uncertainty of the estimated effect for each group.

``` r
ranks <- expectedRank(m1, groupFctr = "d")
head(ranks)
#>      d (Intercept) (Intercept)_var       ER pctER
#> 1 1866   1.2553613     0.012755634 1123.806   100
#> 2 1258   1.1674852     0.034291228 1115.766    99
#> 3  240   1.0933372     0.008761218 1115.090    99
#> 4   79   1.0998653     0.023095979 1112.315    99
#> 5  676   1.0169070     0.026562174 1101.553    98
#> 6   66   0.9568607     0.008602823 1098.049    97
```

Effect Simulation
-----------------

It can still be difficult to interpret the results of LMM and GLMM models, especially the relative influence of varying parameters on the predicted outcome. This is where the `REimpact` and the `wiggle` functions in `merTools` can be handy.

``` r
impSim <- REimpact(m1, InstEval[7, ], groupFctr = "d", breaks = 5, 
                   n.sims = 300, level = 0.9)
impSim
#>   case bin   AvgFit     AvgFitSE nobs
#> 1    1   1 2.799281 3.174322e-04  193
#> 2    1   2 3.278761 6.471135e-05  240
#> 3    1   3 3.573977 5.674664e-05  254
#> 4    1   4 3.853878 5.975449e-05  265
#> 5    1   5 4.244057 1.885801e-04  176
```

The result of `REimpact` shows the change in the `yhat` as the case we supplied to `newdata` is moved from the first to the fifth quintile in terms of the magnitude of the group factor coefficient. We can see here that the individual professor effect has a strong impact on the outcome variable. This can be shown graphically as well:

``` r
library(ggplot2)
ggplot(impSim, aes(x = factor(bin), y = AvgFit, ymin = AvgFit - 1.96*AvgFitSE, 
                   ymax = AvgFit + 1.96*AvgFitSE)) + 
  geom_pointrange() + theme_bw() + labs(x = "Bin of `d` term", y = "Predicted Fit")
```

![](README-reImpactplot-1.png)

Here the standard error is a bit different -- it is the weighted standard error of the mean effect within the bin. It does not take into account the variability within the effects of each observation in the bin -- accounting for this variation will be a future addition to `merTools`.
