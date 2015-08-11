# merTools

[![Travis-CI Build Status](https://travis-ci.org/jknowles/merTools.png?branch=master)](https://travis-ci.org/jknowles/merTools)
[![Coverage Status](https://coveralls.io/repos/jknowles/merTools/badge.svg?branch=master)](https://coveralls.io/r/jknowles/merTools?branch=master)

Working with generalized linear mixed models (GLMM) and linear mixed models (LMM) 
has become increasingly easy with advances in the `lme4` package. 
As we have found ourselves using these models more and more within our work, we, 
the authors, have developed a set of tools for simplifying and speeding up common 
tasks for interacting with `merMod` objects from `lme4`. This package provides 
those tools. 

The easiest way to demo the features of this application is to use the bundled 
Shiny application which launches a number of the metrics here to aide in exploring 
the model. To do this:

```{r}
devtools::install_github("jknowles/merTools")
library(merTools)
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)
shinyMer(m1, simData = InstEval[1:100, ]) # just try the first 100 rows of data
```

On the first tab, the function presents the prediction intervals for the data 
selected by user which are calculated using the `predictInterval` function 
within the package. This function calculates prediction intervals quickly by 
sampling from the simulated distribution of the fixed effect and random effect 
terms and combining these simulated estimates to produce a distribution of 
predictions for each observation. This allows prediction intervals to be generated 
from very large models where the use of `bootMer` would not be feasible 
computationally. 

On the next tab the distribution of the fixed effect and group-level effects 
is depicted on confidence interval plots. These are useful for diagnostics and 
provide a way to inspect the relative magnitudes of various parameters. This 
tab makes use of four related functions in `merTools`: `FEsim`, `plotFEsim`, 
`REsim` and `plotREsim` which are available to be used on their own as well. 

On the third tab are some convenient ways to show the influence or magnitude of 
effects by leveraging the power of `predictInterval`. For each case, up to 12, 
in the selected data type, the user can view the impact of changing either one 
of the fixed effect or one of the grouping level terms. Using the `REimpact` 
function, each case is simulated with the model's prediction if all else was 
held equal, but the observation was moved through the distribution of the 
fixed effect or the random effect term. This is plotted on the scale of the 
dependent variable, which allows the user to compare the magnitude of effects 
across variables, and also between models on the same data. 
