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
merTools
========

A package for getting the most of our multilevel models in R

by Jared E. Knowles and Carl Frederick

Working with generalized linear mixed models (GLMM) and linear mixed
models (LMM) has become increasingly easy with advances in the `lme4`
package. As we have found ourselves using these models more and more
within our work, we, the authors, have developed a set of tools for
simplifying and speeding up common tasks for interacting with `merMod`
objects from `lme4`. This package provides those tools.

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

-   Improve handling of formulas. If the original `merMod` has functions
    specified in the formula, the `draw` and `wiggle` functions will
    check for this and attempt to respect these variable
    transformations. Where this is not possible a warning will be
    issued. Most common transformations are respected as long as the the
    original variable is passed untransformed to the model.
-   Change the calculations of the residual variance. Previously
    residual variance was used to inflate both the variance around the
    fixed parameters and around the predicted values themselves. This
    was incorrect and resulted in overly conservative estimates. Now the
    residual variance is appropriately only used around the final
    predictions
-   New option for `predictInterval` that allows the user to return the
    full interval, the fixed component, the random component, or the
    fixed and each random component separately for each observation
-   Fixed a bug with slope+intercept random terms that caused a
    miscalculation of the random component
-   Add comparison to `rstanarm` to the Vignette
-   Make `expectedRank` output more `tidy` like and allow function to
    calculate expected rank for all terms at once
    -   Note, this breaks the API by changing the names of the columns
        in the output of this function
-   Remove tests that test for timing to avoid issues with R-devel JIT
    compiler
-   Remove `plyr` and replace with `dplyr`
-   Fix issue \#62 `varList` will now throw an error if `==` is used
    instead of `=`
-   Fix issue \#54 `predictInterval` did not included random effects in
    calculations when `newdata` had more than 1000 rows and/or user
    specified `parallel=TRUE`. Note: fix was to disable the `.paropts`
    option for `predictInterval` … user can still specify for
    *temporary* backward compatibility but this should be either removed
    or fixed in the permanent solution.
-   Fix issue \#53 about problems with `predictInterval` when only
    specific levels of a grouping factor are in `newdata` with the colon
    specification of interactions
-   Fix issue \#52 ICC wrong calculations … we just needed to square the
    standard deviations that we pulled

See [NEWS.md](https://github.com/jknowles/merTools/blob/master/NEWS.md)
for more details.

Shiny App and Demo
------------------

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

With `predictInterval` we obtain predictions that are more like the
standard objects produced by `lm` and `glm`:

``` r
#predictInterval(m1, newdata = InstEval[1:10, ]) # all other parameters are optional
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 500, level = 0.9, 
                stat = 'median')
#>         fit      upr      lwr
#> 1  3.217419 5.222041 1.209191
#> 2  3.160320 5.153260 1.077464
#> 3  3.350286 5.334528 1.410776
#> 4  3.167147 5.115386 1.018051
#> 5  3.364984 5.545141 1.344050
#> 6  3.171919 5.269084 1.270331
#> 7  4.242407 5.990232 2.310984
#> 8  3.823723 5.763732 1.860317
#> 9  3.635615 5.644665 1.864689
#> 10 3.208120 5.400739 1.291864
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
#>      effect         fit      upr        lwr obs
#> 1  combined  3.10838502 5.123547  0.9913693   1
#> 2  combined  3.00832927 5.411406  1.0113278   2
#> 3  combined  3.22682784 5.371182  1.2233180   3
#> 4  combined  3.19987923 5.302901  1.0819426   4
#> 5  combined  3.33376041 5.019445  1.3380125   5
#> 6  combined  3.31108813 5.238900  1.2690864   6
#> 7  combined  4.24593858 5.933879  2.1953436   7
#> 8  combined  3.75983116 5.871957  1.8359312   8
#> 9  combined  3.75325587 5.760148  2.0332515   9
#> 10 combined  3.35497451 5.371630  1.4761959  10
#> 11        s  0.34207622 2.567831 -1.6782520   1
#> 12        s  0.19641638 2.360672 -1.7598827   2
#> 13        s -0.03755823 2.236397 -1.4514586   3
#> 14        s  0.07338433 1.931521 -2.0330254   4
#> 15        s -0.13145826 2.035637 -2.4907082   5
#> 16        s -0.06332697 2.058394 -1.9137445   6
#> 17        s  0.39965259 2.207372 -1.8133160   7
#> 18        s  0.35494937 2.363900 -1.5459471   8
#> 19        s  0.18062186 2.171327 -1.7177457   9
#> 20        s  0.39414930 2.437838 -1.4663560  10
#> 21        d -0.19850572 1.868795 -2.0216443   1
#> 22        d -0.23150868 1.927724 -2.0879271   2
#> 23        d  0.12322658 1.973800 -1.9759552   3
#> 24        d -0.21774838 1.878971 -2.2175861   4
#> 25        d  0.13925047 2.121110 -1.9064912   5
#> 26        d  0.07385768 2.105597 -1.9957663   6
#> 27        d  0.54384220 2.373974 -0.9731128   7
#> 28        d  0.33809540 2.342702 -1.8649070   8
#> 29        d  0.32012926 2.260106 -1.6999460   9
#> 30        d -0.29980274 1.521508 -2.2328740  10
#> 31    fixed  3.23207074 5.329527  1.0150544   1
#> 32    fixed  3.11042822 5.304020  1.2880103   2
#> 33    fixed  3.22825199 5.200519  1.2316590   3
#> 34    fixed  3.04045774 4.934492  1.1924062   4
#> 35    fixed  3.20004158 4.977114  1.0328596   5
#> 36    fixed  3.47221763 5.397403  1.3794685   6
#> 37    fixed  3.12964189 5.018796  1.2485740   7
#> 38    fixed  3.44092080 5.523848  1.4783921   8
#> 39    fixed  3.33680051 5.459012  1.2086844   9
#> 40    fixed  3.32120361 5.309080  1.5144517  10
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

![](man/figures/README_unnamed-chunk-8-1.png)

We can also investigate the makeup of the prediction for each
observation.

``` r
ggplot(plotdf[plotdf$obs < 6,], 
       aes(x = effect, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  facet_grid(residVar~obs) + theme_bw()
```

![](man/figures/README_unnamed-chunk-9-1.png)

Plotting
--------

`merTools` also provides functionality for inspecting `merMod` objects
visually. The easiest are getting the posterior distributions of both
fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22261732  3.22132431 0.01781311
#> 2    service1 -0.07299497 -0.07143995 0.01460292
#> 3   lectage.L -0.18689932 -0.18481094 0.01531841
#> 4   lectage.Q  0.02464238  0.02398260 0.01141715
#> 5   lectage.C -0.02454276 -0.02491987 0.01182860
#> 6   lectage^4 -0.01952495 -0.02049488 0.01398923
```

And we can also plot this:

``` r
plotFEsim(FEsim(m1, n.sims = 100), level = 0.9, stat = 'median', intercept = FALSE)
```

![](man/figures/README_FEsimPlot-1.png)

We can also quickly make caterpillar plots for the random-effect terms:

``` r
reSims <- REsim(m1, n.sims = 100)
head(reSims)
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.12289317  0.07537743 0.3265776
#> 2         s       2 (Intercept) -0.06294209 -0.05879910 0.3496404
#> 3         s       3 (Intercept)  0.30066403  0.26274916 0.2963047
#> 4         s       4 (Intercept)  0.23317160  0.26947709 0.2900297
#> 5         s       5 (Intercept)  0.05470018  0.06032355 0.3117877
#> 6         s       6 (Intercept)  0.08529472  0.07354142 0.2396258
```

``` r
plotREsim(REsim(m1, n.sims = 100), stat = 'median', sd = TRUE)
```

![](man/figures/README_reSimplot-1.png)

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

![](man/figures/README_unnamed-chunk-13-1.png)

Effect Simulation
-----------------

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
#> 1    1   1 2.776054 3.007958e-04  193
#> 2    1   2 3.247641 7.012499e-05  240
#> 3    1   3 3.541909 5.622570e-05  254
#> 4    1   4 3.825994 6.723618e-05  265
#> 5    1   5 4.211194 1.849960e-04  176
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

![](man/figures/README_reImpactplot-1.png)

Here the standard error is a bit different – it is the weighted standard
error of the mean effect within the bin. It does not take into account
the variability within the effects of each observation in the bin –
accounting for this variation will be a future addition to `merTools`.

Explore Substantive Impacts
---------------------------

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
variables `situ` and `btype` respectively.

``` r
plotdf <- predictInterval(fmVA, newdata = newData, type = "probability", 
            stat = "median", n.sims = 1000)
plotdf <- cbind(plotdf, newData)

ggplot(plotdf, aes(y = fit, x = Anger, color = btype, group = btype)) + 
  geom_point() + geom_smooth(aes(color = btype), method = "lm") + 
  facet_wrap(~situ) + theme_bw() +
  labs(y = "Predicted Probability")
```

![](man/figures/README_substImpactPredict-1.png)
