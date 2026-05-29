# Analyzing Imputed Data with Multilevel Models and merTools

## Introduction

Multilevel models are valuable in a wide array of problem areas that
involve non-experimental, or observational data. In many of these cases
the data on individual observations may be incomplete. In these
situations, the analyst may turn to one of many methods for filling in
missing data depending on the specific problem at hand, disciplinary
norms, and prior research.

One of the most common cases is to use multiple imputation. Multiple
imputation involves fitting a model to the data and estimating the
missing values for observations. For details on multiple imputation, and
a discussion of some of the main implementations in R, look at the
documentation and vignettes for the `mice` and `Amelia` packages.

The key difficulty multiple imputation creates for users of multilevel
models is that the result of multiple imputation is K replicated
datasets corresponding to different estimated values for the missing
data in the original dataset.

For the purposes of this vignette, I will describe how to use one flavor
of multiple imputation and the function in `merTools` to obtain
estimates from a multilevel model in the presence of missing and
multiply imputed data.

## Missing Data and its Discontents

To demonstrate this workflow, we will use the `hsb` dataset in the
`merTools` package which includes data on the math achievement of a wide
sample of students nested within schools. The data has no missingness,
so first we will simulate some missing data.

``` r

data(hsb)

# Create a function to randomly assign NA values

add_NA <- function(x, prob){
  z <- rbinom(length(x), 1, prob = prob)
  x[z==1] <- NA
  return(x)
}

hsb$minority <- add_NA(hsb$minority, prob = 0.05)
table(is.na(hsb$minority))
#> 
#> FALSE  TRUE 
#>  6832   353

hsb$female <- add_NA(hsb$female, prob = 0.05)
table(is.na(hsb$female))
#> 
#> FALSE  TRUE 
#>  6852   333

hsb$ses <- add_NA(hsb$ses, prob = 0.05)
table(is.na(hsb$ses))
#> 
#> FALSE  TRUE 
#>  6825   360

hsb$size <- add_NA(hsb$size, prob = 0.05)
table(is.na(hsb$size))
#> 
#> FALSE  TRUE 
#>  6815   370
```

``` r

# Load imputation library
library(Amelia)
# Declare the variables to include in the imputation data
varIndex <- names(hsb)
# Declare ID variables to be excluded from imputation
IDS <- c("schid", "meanses")
# Imputate
impute.out <- amelia(hsb[, varIndex], idvars = IDS,
                         noms = c("minority", "female"),
                         m = 5)
#> -- Imputation 1 --
#> 
#>   1  2  3  4
#> 
#> -- Imputation 2 --
#> 
#>   1  2  3
#> 
#> -- Imputation 3 --
#> 
#>   1  2  3
#> 
#> -- Imputation 4 --
#> 
#>   1  2  3
#> 
#> -- Imputation 5 --
#> 
#>   1  2  3
summary(impute.out)
#> 
#> Amelia output with 5 imputed datasets.
#> Return code:  1 
#> Message:  Normal EM convergence. 
#> 
#> Chain Lengths:
#> --------------
#> Imputation 1:  4
#> Imputation 2:  3
#> Imputation 3:  3
#> Imputation 4:  3
#> Imputation 5:  3
#> 
#> Rows after Listwise Deletion:  5870 
#> Rows after Imputation:  7185 
#> Patterns of missingness in the data:  14 
#> 
#> Fraction Missing for original variables: 
#> -----------------------------------------
#> 
#>          Fraction Missing
#> schid          0.00000000
#> minority       0.04913013
#> female         0.04634656
#> ses            0.05010438
#> mathach        0.00000000
#> size           0.05149617
#> schtype        0.00000000
#> meanses        0.00000000
```

``` r

# Amelia is not available so let's just boostrap resample our data
impute.out <- vector(mode = "list", 5)

for (i in 1:5) {
  impute.out[[i]] <- hsb[sample(nrow(hsb), nrow(hsb), replace = TRUE), ]
}

# Declare the variables to include in the imputation data
summary(impute.out)
```

## Fitting and Summarizing a Model List

Fitting a model is very similar

``` r

fmla <- "mathach ~ minority + female + ses + meanses + (1 + ses|schid)"
mod <- lmer(fmla, data = hsb)
if(amelia_eval) {
  modList <- lmerModList(fmla, data = impute.out$imputations)
} else {
  # Use bootstrapped data instead
  modList <- lmerModList(fmla, data = impute.out)
}
```

The resulting object `modList` is a list of `merMod` objects the same
length as the number of imputation datasets. This object is assigned the
class of `merModList` and `merTools` provides some convenience functions
for reporting the results of this object.

Using this, we can directly compare the model fit with missing data
excluded to the aggregate from the imputed models:

``` r

fixef(mod) # model with dropped missing
#> (Intercept)    minority      female         ses     meanses 
#>   14.115316   -2.926850   -1.295694    1.927404    3.111054
fixef(modList)
#> (Intercept)    minority      female         ses     meanses 
#>   14.010651   -2.703166   -1.151025    1.907475    3.129565
```

``` r

VarCorr(mod) # model with dropped missing
#>  Groups   Name        Std.Dev. Corr   
#>  schid    (Intercept) 1.51731         
#>           ses         0.47516  -0.823 
#>  Residual             5.97123
VarCorr(modList) # aggregate of imputed models
#> $stddev
#> $stddev$schid
#> (Intercept)         ses 
#>   1.5146901   0.4775618 
#> 
#> 
#> $correlation
#> $correlation$schid
#>             (Intercept)        ses
#> (Intercept)   1.0000000 -0.8056484
#> ses          -0.8056484  1.0000000
```

If you want to inspect the individual models, or you do not like taking
the mean across the imputation replications, you can take the
`merModList` apart easily:

``` r

lapply(modList, fixef)
#> $imp1
#> (Intercept)    minority      female         ses     meanses 
#>   13.998828   -2.706329   -1.131058    1.919838    3.145374 
#> 
#> $imp2
#> (Intercept)    minority      female         ses     meanses 
#>   14.021477   -2.725340   -1.163246    1.873299    3.164971 
#> 
#> $imp3
#> (Intercept)    minority      female         ses     meanses 
#>   13.960872   -2.509095   -1.156671    1.927857    3.146464 
#> 
#> $imp4
#> (Intercept)    minority      female         ses     meanses 
#>   14.047832   -2.835816   -1.133556    1.916885    3.061399 
#> 
#> $imp5
#> (Intercept)    minority      female         ses     meanses 
#>   14.024244   -2.739249   -1.170593    1.899499    3.129619
```

And, you can always operate on any single element of the list:

``` r

fixef(modList[[1]])
#> (Intercept)    minority      female         ses     meanses 
#>   13.998828   -2.706329   -1.131058    1.919838    3.145374
fixef(modList[[2]])
#> (Intercept)    minority      female         ses     meanses 
#>   14.021477   -2.725340   -1.163246    1.873299    3.164971
```

## Output of a Model List

``` r

print(modList)
#> $imp1
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: mathach ~ minority + female + ses + meanses + (1 + ses | schid)
#>    Data: d
#> 
#> REML criterion at convergence: 46335.2
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2537 -0.7166  0.0341  0.7600  2.9397 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr  
#>  schid    (Intercept)  2.2980  1.5159         
#>           ses          0.2297  0.4793   -0.90 
#>  Residual             35.8645  5.9887         
#> Number of obs: 7185, groups:  schid, 160
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  13.9988     0.1725  81.139
#> minority     -2.7063     0.1978 -13.684
#> female       -1.1311     0.1577  -7.171
#> ses           1.9198     0.1148  16.730
#> meanses       3.1454     0.3493   9.005
#> 
#> Correlation of Fixed Effects:
#>          (Intr) minrty female ses   
#> minority -0.323                     
#> female   -0.480  0.007              
#> ses      -0.276  0.144  0.046       
#> meanses  -0.112  0.122  0.022 -0.247
#> 
#> $imp2
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: mathach ~ minority + female + ses + meanses + (1 + ses | schid)
#>    Data: d
#> 
#> REML criterion at convergence: 46346.6
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2317 -0.7188  0.0317  0.7662  2.9102 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr  
#>  schid    (Intercept)  2.321   1.5235         
#>           ses          0.223   0.4722   -0.88 
#>  Residual             35.915   5.9929         
#> Number of obs: 7185, groups:  schid, 160
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  14.0215     0.1734  80.844
#> minority     -2.7253     0.1975 -13.798
#> female       -1.1632     0.1578  -7.372
#> ses           1.8733     0.1139  16.452
#> meanses       3.1650     0.3507   9.025
#> 
#> Correlation of Fixed Effects:
#>          (Intr) minrty female ses   
#> minority -0.324                     
#> female   -0.485  0.018              
#> ses      -0.266  0.142  0.037       
#> meanses  -0.113  0.124  0.029 -0.241
#> 
#> $imp3
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: mathach ~ minority + female + ses + meanses + (1 + ses | schid)
#>    Data: d
#> 
#> REML criterion at convergence: 46360.8
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2194 -0.7200  0.0327  0.7646  2.9237 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr  
#>  schid    (Intercept)  2.234   1.4947         
#>           ses          0.275   0.5244   -0.56 
#>  Residual             35.956   5.9963         
#> Number of obs: 7185, groups:  schid, 160
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  13.9609     0.1719  81.198
#> minority     -2.5091     0.1978 -12.688
#> female       -1.1567     0.1585  -7.299
#> ses           1.9279     0.1161  16.611
#> meanses       3.1465     0.3538   8.894
#> 
#> Correlation of Fixed Effects:
#>          (Intr) minrty female ses   
#> minority -0.327                     
#> female   -0.488  0.015              
#> ses      -0.208  0.144  0.045       
#> meanses  -0.096  0.127  0.025 -0.242
#> 
#> $imp4
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: mathach ~ minority + female + ses + meanses + (1 + ses | schid)
#>    Data: d
#> 
#> REML criterion at convergence: 46313.8
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2487 -0.7150  0.0369  0.7606  2.9364 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr  
#>  schid    (Intercept)  2.3205  1.5233         
#>           ses          0.1788  0.4229   -0.89 
#>  Residual             35.7569  5.9797         
#> Number of obs: 7185, groups:  schid, 160
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  14.0478     0.1728  81.300
#> minority     -2.8358     0.1969 -14.405
#> female       -1.1336     0.1578  -7.186
#> ses           1.9169     0.1133  16.912
#> meanses       3.0614     0.3529   8.674
#> 
#> Correlation of Fixed Effects:
#>          (Intr) minrty female ses   
#> minority -0.322                     
#> female   -0.480  0.007              
#> ses      -0.249  0.149  0.041       
#> meanses  -0.108  0.121  0.026 -0.248
#> 
#> $imp5
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: mathach ~ minority + female + ses + meanses + (1 + ses | schid)
#>    Data: d
#> 
#> REML criterion at convergence: 46332.3
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2359 -0.7195  0.0323  0.7624  2.9225 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr  
#>  schid    (Intercept)  2.2982  1.516          
#>           ses          0.2391  0.489    -0.80 
#>  Residual             35.8344  5.986          
#> Number of obs: 7185, groups:  schid, 160
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  14.0242     0.1723  81.399
#> minority     -2.7392     0.1973 -13.884
#> female       -1.1706     0.1576  -7.425
#> ses           1.8995     0.1148  16.539
#> meanses       3.1296     0.3513   8.910
#> 
#> Correlation of Fixed Effects:
#>          (Intr) minrty female ses   
#> minority -0.320                     
#> female   -0.479  0.003              
#> ses      -0.257  0.143  0.042       
#> meanses  -0.106  0.122  0.023 -0.243
```

``` r

summary(modList)
#> [1] "Linear mixed model fit by REML"
#> Model family: 
#> lmer(formula = mathach ~ minority + female + ses + meanses + 
#>     (1 + ses | schid), data = d)
#> 
#> Fixed Effects:
#>             estimate std.error statistic          df
#> (Intercept)   14.011     0.174    80.693  113115.704
#> female        -1.151     0.158    -7.277 1080751.936
#> meanses        3.130     0.353     8.862  209086.287
#> minority      -2.703     0.211   -12.803     950.730
#> ses            1.907     0.115    16.583  257649.852
#> 
#> Random Effects:
#> 
#> Error Term Standard Deviations by Level:
#> 
#> schid
#> (Intercept)         ses 
#>       1.515       0.478 
#> 
#> 
#> Error Term Correlations:
#> 
#> schid
#>             (Intercept) ses   
#> (Intercept)  1.000      -0.806
#> ses         -0.806       1.000
#> 
#> 
#> Residual Error = 5.989 
#> 
#> ---Groups
#> number of obs: 7185, groups: schid, 160
#> 
#> Model Fit Stats
#> AIC = 46355.7
#> Residual standard deviation = 5.989
```

``` r

fastdisp(modList)
#> lmer(formula = mathach ~ minority + female + ses + meanses + 
#>     (1 + ses | schid), data = d)
#>             estimate std.error
#> (Intercept)    14.01      0.17
#> female         -1.15      0.16
#> meanses         3.13      0.35
#> minority       -2.70      0.21
#> ses             1.91      0.12
#> 
#> Error terms:
#>  Groups   Name        Std.Dev. Corr  
#>  schid    (Intercept) 1.51           
#>           ses         0.48     -0.90 
#>  Residual             5.99           
#> ---
#> number of obs: 7185, groups: schid, 160
#> AIC = 46355.7---
```

The standard errors reported for the model list include a correction,
Rubin’s correction (see documentation), which adjusts for the within and
between imputation set variance as well.

## Specific Model Information Summaries

``` r

modelRandEffStats(modList)
#>                   term    group   estimate   std.error
#> 1 cor__(Intercept).ses    schid -0.8056484 0.141646398
#> 2      sd__(Intercept)    schid  1.5146901 0.011792677
#> 3      sd__Observation Residual  5.9887748 0.006394734
#> 4              sd__ses    schid  0.4775618 0.036568245
modelFixedEff(modList)
#>          term  estimate std.error  statistic           df
#> 1 (Intercept) 14.010651 0.1736290  80.693025  113115.7037
#> 2      female -1.151025 0.1581831  -7.276536 1080751.9364
#> 3     meanses  3.129565 0.3531353   8.862227  209086.2872
#> 4    minority -2.703166 0.2111383 -12.802819     950.7301
#> 5         ses  1.907475 0.1150281  16.582694  257649.8518
VarCorr(modList)
#> $stddev
#> $stddev$schid
#> (Intercept)         ses 
#>   1.5146901   0.4775618 
#> 
#> 
#> $correlation
#> $correlation$schid
#>             (Intercept)        ses
#> (Intercept)   1.0000000 -0.8056484
#> ses          -0.8056484  1.0000000
```

### Diagnostics of List Components

``` r

modelInfo(mod)
#>   n.obs n.lvls      AIC    sigma
#> 1  6177      1 39832.18 5.971227
```

Let’s apply this to our model list.

``` r

lapply(modList, modelInfo)
#> $imp1
#>   n.obs n.lvls      AIC    sigma
#> 1  7185      1 46353.21 5.988695
#> 
#> $imp2
#>   n.obs n.lvls     AIC    sigma
#> 1  7185      1 46364.6 5.992938
#> 
#> $imp3
#>   n.obs n.lvls     AIC    sigma
#> 1  7185      1 46378.8 5.996347
#> 
#> $imp4
#>   n.obs n.lvls      AIC   sigma
#> 1  7185      1 46331.78 5.97971
#> 
#> $imp5
#>   n.obs n.lvls      AIC    sigma
#> 1  7185      1 46350.26 5.986183
```

### Model List Generics

``` r

summary(modList)
#> [1] "Linear mixed model fit by REML"
#> Model family: 
#> lmer(formula = mathach ~ minority + female + ses + meanses + 
#>     (1 + ses | schid), data = d)
#> 
#> Fixed Effects:
#>             estimate std.error statistic          df
#> (Intercept)   14.011     0.174    80.693  113115.704
#> female        -1.151     0.158    -7.277 1080751.936
#> meanses        3.130     0.353     8.862  209086.287
#> minority      -2.703     0.211   -12.803     950.730
#> ses            1.907     0.115    16.583  257649.852
#> 
#> Random Effects:
#> 
#> Error Term Standard Deviations by Level:
#> 
#> schid
#> (Intercept)         ses 
#>       1.515       0.478 
#> 
#> 
#> Error Term Correlations:
#> 
#> schid
#>             (Intercept) ses   
#> (Intercept)  1.000      -0.806
#> ses         -0.806       1.000
#> 
#> 
#> Residual Error = 5.989 
#> 
#> ---Groups
#> number of obs: 7185, groups: schid, 160
#> 
#> Model Fit Stats
#> AIC = 46355.7
#> Residual standard deviation = 5.989
```

``` r

modelFixedEff(modList)
#>          term  estimate std.error  statistic           df
#> 1 (Intercept) 14.010651 0.1736290  80.693025  113115.7037
#> 2      female -1.151025 0.1581831  -7.276536 1080751.9364
#> 3     meanses  3.129565 0.3531353   8.862227  209086.2872
#> 4    minority -2.703166 0.2111383 -12.802819     950.7301
#> 5         ses  1.907475 0.1150281  16.582694  257649.8518
```

``` r

ranef(modList)
#> $schid
#>       (Intercept)           ses
#> 1224 -0.169934903  0.0484630054
#> 1288 -0.007840719  0.0020184342
#> 1296 -0.120121384  0.0316003121
#> 1308  0.110260016 -0.0345814029
#> 1317  0.068892627 -0.0236912696
#> 1358 -0.296774372  0.0904682545
#> 1374 -0.405072947  0.1220606747
#> 1433  0.318400202 -0.0830043713
#> 1436  0.259394301 -0.0698354596
#> 1461 -0.078477863  0.0391032731
#> 1462  0.343550954 -0.1048708761
#> 1477  0.035716104 -0.0164384112
#> 1499 -0.317868388  0.0968881791
#> 1637 -0.090390950  0.0269646899
#> 1906  0.055768347 -0.0169288276
#> 1909 -0.065472718  0.0169395570
#> 1942  0.196755113 -0.0558813645
#> 1946 -0.095505668  0.0402309899
#> 2030 -0.402888189  0.1036469112
#> 2208 -0.025120302  0.0114460657
#> 2277  0.337550172 -0.1136157208
#> 2305  0.582441616 -0.1854428413
#> 2336  0.134604118 -0.0366408100
#> 2458  0.241996781 -0.0646483272
#> 2467 -0.226680562  0.0625779063
#> 2526  0.439556485 -0.1276085840
#> 2626  0.012854549  0.0012862800
#> 2629  0.360692714 -0.1075798863
#> 2639  0.120441776 -0.0423993743
#> 2651 -0.406708936  0.1242891680
#> 2655  0.615656442 -0.1693982670
#> 2658 -0.261375232  0.0733533865
#> 2755  0.111844309 -0.0347401814
#> 2768 -0.287438525  0.0866945732
#> 2771  0.007098193  0.0055279386
#> 2818 -0.017209992  0.0099641840
#> 2917  0.170991702 -0.0565713328
#> 2990  0.482050460 -0.1316049132
#> 2995 -0.222728540  0.0531943679
#> 3013 -0.126917295  0.0408401321
#> 3020  0.054081114 -0.0155796367
#> 3039  0.244154050 -0.0654082046
#> 3088 -0.017186488 -0.0004443775
#> 3152 -0.054853716  0.0211548556
#> 3332 -0.257656248  0.0698106244
#> 3351 -0.420868035  0.1138012619
#> 3377  0.174422620 -0.0631883751
#> 3427  0.875949018 -0.2551272985
#> 3498 -0.009415170 -0.0069280841
#> 3499 -0.138316005  0.0318361010
#> 3533 -0.145407439  0.0316613782
#> 3610  0.313525367 -0.0790018448
#> 3657 -0.086102109  0.0280277586
#> 3688 -0.030130748  0.0074381151
#> 3705 -0.385929782  0.1023309955
#> 3716  0.035108793  0.0052430806
#> 3838  0.479257816 -0.1411916170
#> 3881 -0.296850514  0.0823500523
#> 3967 -0.045346466  0.0181219299
#> 3992  0.099519612 -0.0356231158
#> 3999 -0.052640860  0.0216319988
#> 4042 -0.197893269  0.0488311886
#> 4173 -0.088350067  0.0301256054
#> 4223  0.266744136 -0.0778116717
#> 4253  0.031749842 -0.0255555401
#> 4292  0.521310978 -0.1594846494
#> 4325  0.009310789  0.0050361203
#> 4350 -0.297295651  0.0908865951
#> 4383 -0.282542432  0.0857098301
#> 4410 -0.074926191  0.0232212238
#> 4420  0.184666268 -0.0496780305
#> 4458 -0.021676941  0.0017713352
#> 4511  0.217381792 -0.0672136642
#> 4523 -0.259948575  0.0760573128
#> 4530  0.017403274 -0.0020254277
#> 4642  0.125767623 -0.0273296529
#> 4868 -0.236776899  0.0564516730
#> 4931 -0.164323177  0.0369784096
#> 5192 -0.214861224  0.0573871806
#> 5404 -0.208359199  0.0449287006
#> 5619 -0.054056780  0.0318423636
#> 5640  0.064664160 -0.0126144003
#> 5650  0.489479118 -0.1407784263
#> 5667 -0.283509096  0.0857569790
#> 5720  0.097236634 -0.0245184608
#> 5761  0.119344934 -0.0242388711
#> 5762 -0.071197537  0.0162747328
#> 5783 -0.081933178  0.0247260166
#> 5815 -0.187340951  0.0555854868
#> 5819 -0.315919366  0.0849850257
#> 5838 -0.029283350  0.0077877641
#> 5937  0.050432229 -0.0176749438
#> 6074  0.360234770 -0.1029339227
#> 6089  0.224408433 -0.0604356535
#> 6144 -0.291204156  0.0851950899
#> 6170  0.290597283 -0.0742598668
#> 6291  0.205954836 -0.0561388343
#> 6366  0.202510633 -0.0594848229
#> 6397  0.168407017 -0.0444952531
#> 6415 -0.066353091  0.0252905737
#> 6443 -0.074735917  0.0102799339
#> 6464 -0.026313523  0.0051057369
#> 6469  0.303898400 -0.0785744370
#> 6484  0.114243188 -0.0426948964
#> 6578  0.312923327 -0.0891869407
#> 6600 -0.245760181  0.0803257940
#> 6808 -0.303255196  0.0811582973
#> 6816  0.180927859 -0.0524687532
#> 6897  0.034103208 -0.0006581543
#> 6990 -0.285854793  0.0721911619
#> 7011  0.031240836 -0.0027379460
#> 7101 -0.121583610  0.0331155054
#> 7172 -0.185365390  0.0481045833
#> 7232 -0.046127481  0.0263432687
#> 7276 -0.079475511  0.0312825132
#> 7332  0.055182822 -0.0158540995
#> 7341 -0.284399889  0.0720398869
#> 7342  0.064212989 -0.0189717669
#> 7345 -0.252890323  0.0805535321
#> 7364  0.290202961 -0.0850729187
#> 7635  0.073112761 -0.0172150034
#> 7688  0.619349977 -0.1792330219
#> 7697  0.095639141 -0.0238826906
#> 7734  0.014643894  0.0086617659
#> 7890 -0.275256456  0.0715160781
#> 7919 -0.164444669  0.0508069476
#> 8009 -0.225761029  0.0618499297
#> 8150  0.075671140 -0.0279229621
#> 8165  0.171725217 -0.0491859344
#> 8175  0.123498307 -0.0385463505
#> 8188 -0.150869785  0.0448245563
#> 8193  0.514473612 -0.1438987815
#> 8202 -0.239384601  0.0729096036
#> 8357  0.197473961 -0.0562790086
#> 8367 -0.776217194  0.2141058702
#> 8477  0.036517616 -0.0007453225
#> 8531 -0.219928087  0.0626051566
#> 8627 -0.352196603  0.0942615438
#> 8628  0.595899916 -0.1704784683
#> 8707 -0.076777352  0.0261532425
#> 8775 -0.190797237  0.0470764418
#> 8800  0.007145467 -0.0025477033
#> 8854 -0.575354318  0.1649307995
#> 8857  0.269589997 -0.0787290632
#> 8874  0.154033890 -0.0360728209
#> 8946 -0.158654369  0.0475701902
#> 8983 -0.124362965  0.0295975945
#> 9021 -0.261325833  0.0705960802
#> 9104  0.004559101 -0.0033644879
#> 9158 -0.292540123  0.0920067855
#> 9198  0.347708797 -0.0975904342
#> 9225  0.018927189 -0.0005751321
#> 9292  0.274423361 -0.0821335974
#> 9340 -0.022945413  0.0092996609
#> 9347 -0.062717973  0.0244423946
#> 9359  0.012302584 -0.0180649618
#> 9397 -0.520016492  0.1450591625
#> 9508  0.085206319 -0.0193617510
#> 9550 -0.267133581  0.0801879186
#> 9586 -0.153549872  0.0419366968
```

## Cautions and Notes

Often it is desirable to include aggregate values in the level two or
level three part of the model such as level 1 SES and level 2 mean SES
for the group. In cases where there is missingness in either the level 1
SES values, or in the level 2 mean SES values, caution and careful
thought need to be given to how to proceed with the imputation routine.
