pkgname <- "merTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('merTools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("FEsim")
### * FEsim

flush(stderr()); flush(stdout())

### Name: FEsim
### Title: Simulate fixed effects from merMod 'FEsim' simulates fixed
###   effects from merMod object posterior distributions
### Aliases: FEsim

### ** Examples

require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fe2 <- FEsim(m2, 25)
head(fe2)



cleanEx()
nameEx("ICC")
### * ICC

flush(stderr()); flush(stdout())

### Name: ICC
### Title: Calculate the intraclass correlation using mixed effect models
### Aliases: ICC

### ** Examples

data(sleepstudy)
ICC(outcome = "Reaction", group = "Subject", data = sleepstudy)



cleanEx()
nameEx("REcorrExtract")
### * REcorrExtract

flush(stderr()); flush(stdout())

### Name: REcorrExtract
### Title: Extract the correlations between the slopes and the intercepts
###   from a model
### Aliases: REcorrExtract

### ** Examples

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REcorrExtract(fm1)



cleanEx()
nameEx("REextract")
### * REextract

flush(stderr()); flush(stdout())

### Name: REextract
### Title: Extracts random effects
### Aliases: REextract

### ** Examples

require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
rfx <- REextract(m2)
#Note the column names
head(rfx)



cleanEx()
nameEx("REimpact")
### * REimpact

flush(stderr()); flush(stdout())

### Name: REimpact
### Title: Calculate the weighted mean of fitted values for various levels
###   of random effect terms.
### Aliases: REimpact

### ** Examples

#For a one-level random intercept model
require(lme4)
m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
(m1.er <- REimpact(m1, newdata = sleepstudy[1, ], breaks = 2))

#For a one-level random intercept model with multiple random terms
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#ranked by the random slope on Days
(m2.er1 <- REimpact(m2,  newdata = sleepstudy[1, ],
           groupFctr = "Subject", term="Days"))
#ranked by the random intercept
(m2.er2 <- REimpact(m2, newdata = sleepstudy[1, ],
             groupFctr = "Subject", term="int"))
## Not run: 
##D # You can also pass additional arguments to predictInterval through REimpact
##D g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
##D zed <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", n.sims = 50,
##D                 include.resid.var = TRUE)
##D zed2 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "s", n.sims = 50,
##D                  include.resid.var = TRUE)
##D zed3 <- REimpact(g1, newdata = InstEval[9:12, ], groupFctr = "d", breaks = 5,
## End(Not run)




cleanEx()
nameEx("REquantile")
### * REquantile

flush(stderr()); flush(stdout())

### Name: REquantile
### Title: Identify group level associated with RE quantile
### Aliases: REquantile

### ** Examples

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REquantile(fm1, quantile = 0.25, groupFctr = "Subject")
REquantile(fm1, quantile = 0.25, groupFctr = "Subject", term = "Days")



cleanEx()
nameEx("REsdExtract")
### * REsdExtract

flush(stderr()); flush(stdout())

### Name: REsdExtract
### Title: Extract the standard deviation of the random effects from a
###   merMod object
### Aliases: REsdExtract

### ** Examples

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
REsdExtract(fm1)



cleanEx()
nameEx("REsim")
### * REsim

flush(stderr()); flush(stdout())

### Name: REsim
### Title: Simulate random effects from merMod 'REsim' simulates random
###   effects from merMod object posterior distributions
### Aliases: REsim

### ** Examples

require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
re2 <- REsim(m2, 25)
head(re2)



cleanEx()
nameEx("RMSE.merMod")
### * RMSE.merMod

flush(stderr()); flush(stdout())

### Name: RMSE.merMod
### Title: Estimate the Root Mean Squared Error (RMSE) for a lmerMod
### Aliases: RMSE.merMod

### ** Examples

require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
RMSE.merMod(m2)



cleanEx()
nameEx("VarCorr.merModList")
### * VarCorr.merModList

flush(stderr()); flush(stdout())

### Name: VarCorr.merModList
### Title: Extract the variances and correlations for random effects from a
###   merMod list
### Aliases: VarCorr.merModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
VarCorr(mod)



cleanEx()
nameEx("draw")
### * draw

flush(stderr()); flush(stdout())

### Name: draw
### Title: Draw a single observation out of an object matching some
###   criteria
### Aliases: draw draw.merMod

### ** Examples

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# Random case
draw(fm1, type = "random")
# Average
draw(fm1, type = "average")
# Subset
draw(fm1, type = "average", varList = list("Subject" = "308"))




cleanEx()
nameEx("expectedRank")
### * expectedRank

flush(stderr()); flush(stdout())

### Name: expectedRank
### Title: Calculate the expected rank of random coefficients that account
###   for uncertainty.
### Aliases: expectedRank

### ** Examples

#For a one-level random intercept model
require(lme4)
m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
(m1.er <- expectedRank(m1))

#For a one-level random intercept model with multiple random terms
require(lme4)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#ranked by the random slope on Days
(m2.er1 <- expectedRank(m2, term="Days"))
#ranked by the random intercept
(m2.er2 <- expectedRank(m2, term="int"))

## Not run: 
##D #For a two-level model with random intercepts
##D require(lme4)
##D m3 <- lmer(y ~ service * dept + (1|s) + (1|d), InstEval)
##D #Ranked by the random intercept on 's'
##D (m3.er1 <- expectedRank(m3, groupFctr="s", term="Intercept"))
## End(Not run)



cleanEx()
nameEx("fastdisp")
### * fastdisp

flush(stderr()); flush(stdout())

### Name: fastdisp
### Title: fastdisp: faster display of model summaries
### Aliases: fastdisp fastdisp.merMod fastdisp.merModList

### ** Examples

## Not run: 
##D #Compare the time for displaying this modest model
##D require(arm)
##D m1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
##D system.time(display(m1))
##D system.time(fastdisp(m1))
## End(Not run)



cleanEx()
nameEx("fixef.merModList")
### * fixef.merModList

flush(stderr()); flush(stdout())

### Name: fixef.merModList
### Title: Extract fixed-effects estimates for a merModList
### Aliases: fixef.merModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
fixef(mod)



cleanEx()
nameEx("hsb")
### * hsb

flush(stderr()); flush(stdout())

### Name: hsb
### Title: A subset of data from the 1982 High School and Beyond survey
###   used as examples for HLM software
### Aliases: hsb
### Keywords: datasets

### ** Examples

data(hsb)
head(hsb)



cleanEx()
nameEx("merModList")
### * merModList

flush(stderr()); flush(stdout())

### Name: lmerModList
### Title: Apply a multilevel model to a list of data frames
### Aliases: lmerModList blmerModList glmerModList bglmerModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
summary(mod)



cleanEx()
nameEx("modelFixedEff")
### * modelFixedEff

flush(stderr()); flush(stdout())

### Name: modelFixedEff
### Title: Extract averaged fixed effect parameters across a list of merMod
###   objects
### Aliases: modelFixedEff

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelFixedEff(mod)



cleanEx()
nameEx("modelInfo")
### * modelInfo

flush(stderr()); flush(stdout())

### Name: modelInfo
### Title: Extract model information from a merMod
### Aliases: modelInfo

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelInfo(mod[[1]])
lapply(mod, modelInfo)



cleanEx()
nameEx("modelRandEffStats")
### * modelRandEffStats

flush(stderr()); flush(stdout())

### Name: modelRandEffStats
### Title: Extract data.frame of random effect statistics from merMod List
### Aliases: modelRandEffStats

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
modelRandEffStats(mod)



cleanEx()
nameEx("plotFEsim")
### * plotFEsim

flush(stderr()); flush(stdout())

### Name: plotFEsim
### Title: Plot the results of a simulation of the fixed effects
### Aliases: plotFEsim

### ** Examples

 fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
 (p1 <- plotFEsim(FEsim(fm1)))



cleanEx()
nameEx("plotREsim")
### * plotREsim

flush(stderr()); flush(stdout())

### Name: plotREsim
### Title: Plot the results of a simulation of the random effects
### Aliases: plotREsim

### ** Examples

 fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
 (p1 <- plotREsim(REsim(fm1)))
 #Plot just the random effects for the Days slope
 (p2 <- plotREsim(REsim(fm1), facet= list(groupFctr= "Subject", term= "Days")))



cleanEx()
nameEx("predictInterval")
### * predictInterval

flush(stderr()); flush(stdout())

### Name: predictInterval
### Title: Predict from merMod objects with a prediction interval
### Aliases: predictInterval

### ** Examples

m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
regFit <- predict(m1, newdata = sleepstudy[11, ]) # a single value is returned
intFit <- predictInterval(m1, newdata = sleepstudy[11, ]) # bounded values
# Can do glmer
d1 <- cbpp
d1$y <- d1$incidence / d1$size
 gm2 <- glmer(y ~ period + (1 | herd), family = binomial, data = d1,
               nAGQ = 9, weights = d1$size)
 regFit <- predict(gm2, newdata = d1[1:10, ])
 # get probabilities
 regFit <- predict(gm2, newdata = d1[1:10, ], type = "response")
 intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "probability")
 intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "linear.prediction")



cleanEx()
nameEx("print.merModList")
### * print.merModList

flush(stderr()); flush(stdout())

### Name: print.merModList
### Title: Print the results of a merMod list
### Aliases: print.merModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
print(mod)



cleanEx()
nameEx("ranef.merModList")
### * ranef.merModList

flush(stderr()); flush(stdout())

### Name: ranef.merModList
### Title: Extract random-effects estimates for a merModList
### Aliases: ranef.merModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
ranef(mod)



cleanEx()
nameEx("subBoot")
### * subBoot

flush(stderr()); flush(stdout())

### Name: subBoot
### Title: Bootstrap a subset of an lme4 model
### Aliases: subBoot

### ** Examples

## Not run: 
##D (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
##D resultMatrix <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20)
## End(Not run)



cleanEx()
nameEx("summary.merModList")
### * summary.merModList

flush(stderr()); flush(stdout())

### Name: summary.merModList
### Title: Summarize a merMod list
### Aliases: summary.merModList

### ** Examples

sim_list <- replicate(n = 10,
        expr = sleepstudy[sample(row.names(sleepstudy), 180),],
        simplify=FALSE)
fml <- "Reaction ~ Days + (Days | Subject)"
mod <- lmerModList(fml, data = sim_list)
summary(mod)



cleanEx()
nameEx("superFactor")
### * superFactor

flush(stderr()); flush(stdout())

### Name: superFactor
### Title: Create a factor with unobserved levels
### Aliases: superFactor

### ** Examples

regularFactor <- c("A", "B", "C")
regularFactor <- factor(regularFactor)
levels(regularFactor)
# Now make it super
newLevs <- c("D", "E", "F")
regularFactor <- superFactor(regularFactor, fullLev = newLevs)
levels(regularFactor) # now super



cleanEx()
nameEx("thetaExtract")
### * thetaExtract

flush(stderr()); flush(stdout())

### Name: thetaExtract
### Title: Extract theta parameters from a merMod model
### Aliases: thetaExtract

### ** Examples

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
thetaExtract(fm1) #(a numeric vector of the covariance parameters)



cleanEx()
nameEx("wiggle")
### * wiggle

flush(stderr()); flush(stdout())

### Name: wiggle
### Title: Assign an observation to different values
### Aliases: wiggle

### ** Examples

data(iris)
wiggle(iris[3,], varlist = "Sepal.Width", valueslist = list(c(1, 2, 3, 5)))
wiggle(iris[3:5,], "Sepal.Width", valueslist = list(c(1, 2, 3, 5)))
wiggle(iris[3,], c("Sepal.Width", "Petal.Length"), list(c(1,2,3,5), c(3,5,6)))
wiggle(iris[3:5,], c("Sepal.Width", "Petal.Length"), list(c(1,2,3,5), c(3,5,6)))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
