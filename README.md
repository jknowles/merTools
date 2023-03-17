[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/merTools)](https://cran.r-project.org/package=merTools)
[![Downloads](http://cranlogs.r-pkg.org/badges/merTools)](https://cran.r-project.org/package=merTools)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/merTools)](https://cran.r-project.org/package=merTools)

<!-- README.md is generated from README.Rmd. Please edit that file -->

    ## Loading required package: arm

    ## Loading required package: MASS

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.2.2

    ## Loading required package: lme4

    ## Warning: package 'lme4' was built under R version 4.2.2

    ## 
    ## arm (Version 1.13-1, built: 2022-8-25)

    ## Working directory is C:/Users/Jared/GitHub/merTools

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

## merTools 0.6.1 (Spring 2023)

- Maintenance release to keep package listed on CRAN
- Fix a small bug where parallel code path is run twice (#126)
- Update plotting functions to avoid deprecated `aes_string()` calls
  (#127)
- Fix (#115) in description
- Speed up PI using @bbolker pull request (#120)
- Updated package maintainer contact information

### merTools 0.5.0

#### New Features

- `subBoot` now works with `glmerMod` objects as well
- `reMargins` a new function that allows the user to marginalize the
  prediction over breaks in the distribution of random effect
  distributions, see `?reMargins` and the new `reMargins` vignette
  (closes \#73)

#### Bug fixes

- Fixed an issue where known convergence errors were issuing warnings
  and causing the test suite to not work
- Fixed an issue where models with a random slope, no intercept, and no
  fixed term were unable to be predicted (#101)
- Fixed an issue with shinyMer not working with substantive fixed
  effects (#93)

### merTools 0.4.1

#### New Features

- Standard errors reported by `merModList` functions now apply the Rubin
  correction for multiple imputation

#### Bug fixes

- Contribution by Alex Whitworth (@alexWhitworth) adding error checking
  to plotting functions

## Shiny App and Demo

The easiest way to demo the features of this application is to use the
bundled Shiny application which launches a number of the metrics here to
aide in exploring the model. To do this:

    library(merTools)
    m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
    shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data

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
predictInterval(m1, newdata = InstEval[1:10, ], n.sims = 500, level = 0.9, 
                stat = 'median')
#>         fit      upr      lwr
#> 1  3.092041 5.059149 1.102673
#> 2  3.204839 5.194856 1.106432
#> 3  3.416692 5.582825 1.474051
#> 4  3.126408 5.012763 1.155690
#> 5  3.378350 5.432257 1.316291
#> 6  3.189796 5.097584 1.003082
#> 7  4.009440 6.265815 2.062227
#> 8  3.871985 5.665837 1.957945
#> 9  3.786645 5.827380 1.734930
#> 10 3.329288 5.399349 1.118200
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
#>      effect          fit      upr       lwr obs
#> 1  combined  3.195900915 5.052559  1.250406   1
#> 2  combined  3.169189347 5.207359  1.591827   2
#> 3  combined  3.262542635 5.570187  1.459483   3
#> 4  combined  3.115651728 5.084291  1.106097   4
#> 5  combined  3.237728407 5.220338  1.588008   5
#> 6  combined  3.302786120 5.439637  1.423937   6
#> 7  combined  4.189503434 6.185946  2.253511   7
#> 8  combined  3.988773316 5.613463  1.994584   8
#> 9  combined  3.758915193 5.861685  1.803003   9
#> 10 combined  3.304601000 5.458847  1.308125  10
#> 11        s  0.134381011 2.167280 -1.586050   1
#> 12        s -0.001042665 2.467795 -1.915679   2
#> 13        s  0.187569316 2.096917 -1.734940   3
#> 14        s  0.218107920 1.952775 -1.713542   4
#> 15        s  0.025583268 1.720703 -2.101509   5
#> 16        s -0.239272443 1.810596 -2.289241   6
#> 17        s  0.240631453 2.295807 -1.573429   7
#> 18        s  0.412103710 2.369665 -1.412571   8
#> 19        s  0.384798169 2.286189 -1.641814   9
#> 20        s  0.318774764 2.414534 -1.336843  10
#> 21        d -0.326676611 1.562441 -2.358736   1
#> 22        d -0.294823958 1.306738 -2.081158   2
#> 23        d  0.131937169 1.974962 -1.861465   3
#> 24        d -0.266664972 1.706584 -2.028685   4
#> 25        d  0.109290898 1.968827 -1.917052   5
#> 26        d  0.023798126 1.762352 -2.065934   6
#> 27        d  0.449466941 2.182760 -1.398073   7
#> 28        d  0.262878220 2.134833 -1.659361   8
#> 29        d  0.001006160 2.316176 -1.920622   9
#> 30        d -0.364587488 1.394051 -2.558685  10
#> 31    fixed  3.013326701 5.139050  1.214449   1
#> 32    fixed  3.171549982 5.028217  1.412183   2
#> 33    fixed  3.135024555 5.227997  1.378779   3
#> 34    fixed  3.222648677 5.168780  1.177561   4
#> 35    fixed  3.389946040 4.970907  1.500747   5
#> 36    fixed  3.309274692 5.020056  1.572196   6
#> 37    fixed  3.330341791 5.237401  1.627099   7
#> 38    fixed  3.134645303 5.281317  1.274267   8
#> 39    fixed  3.208466137 5.106115  1.198307   9
#> 40    fixed  3.303535683 5.336764  1.361588  10
```

This can lead to some useful plotting:

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.2.2
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
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
```

![](man/figures/README_unnamed-chunk-7-1.png)<!-- -->

We can also investigate the makeup of the prediction for each
observation.

``` r
ggplot(plotdf[plotdf$obs < 6,], 
       aes(x = effect, y = fit, ymin = lwr, ymax = upr)) + 
  geom_pointrange() +
  geom_hline(yintercept = 0, color = I("red"), size = 1.1) +
  facet_grid(residVar~obs) + theme_bw()
```

![](man/figures/README_unnamed-chunk-8-1.png)<!-- -->

## Plotting

`merTools` also provides functionality for inspecting `merMod` objects
visually. The easiest are getting the posterior distributions of both
fixed and random effect parameters.

``` r
feSims <- FEsim(m1, n.sims = 100)
head(feSims)
#>          term        mean      median         sd
#> 1 (Intercept)  3.22414807  3.22415334 0.01799582
#> 2    service1 -0.07144090 -0.07156045 0.01302066
#> 3   lectage.L -0.18562440 -0.18603861 0.01699212
#> 4   lectage.Q  0.02591466  0.02678169 0.01202715
#> 5   lectage.C -0.02571353 -0.02598985 0.01255748
#> 6   lectage^4 -0.02198558 -0.02232974 0.01473034
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
#>   groupFctr groupID        term        mean      median        sd
#> 1         s       1 (Intercept)  0.18980845  0.19710994 0.3255318
#> 2         s       2 (Intercept) -0.02890981 -0.06606199 0.3324246
#> 3         s       3 (Intercept)  0.32110508  0.32929167 0.3062943
#> 4         s       4 (Intercept)  0.26867844  0.24816535 0.3147262
#> 5         s       5 (Intercept)  0.04462127  0.04392913 0.3144734
#> 6         s       6 (Intercept)  0.06556643  0.06913915 0.2198916
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

![](man/figures/README_unnamed-chunk-12-1.png)<!-- -->

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
#> 1    1   1 2.798648 2.830908e-04  193
#> 2    1   2 3.265509 6.090223e-05  240
#> 3    1   3 3.559953 5.726929e-05  254
#> 4    1   4 3.848517 5.614353e-05  265
#> 5    1   5 4.237884 2.037161e-04  176
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
#> `geom_smooth()` using formula = 'y ~ x'
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

![](man/figures/README_unnamed-chunk-14-1.png)<!-- -->
