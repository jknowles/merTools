# -----------------------------------------------------
# Test framework includes tests for multiple intercepts and
# multiple slopes to ensure extraction of random effects
# works in these scenarios
#-------------------------------------------------------

set.seed(51315)
library(lme4)
data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
data(grouseticks_agg)
grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")

# Build out models
form <- TICKS ~ YEAR + HEIGHT +(1|BROOD) + (1|INDEX) + (1|LOCATION)
glmer3Lev  <- glmer(form, family="poisson",data=grouseticks,
                     control = glmerControl(optimizer="Nelder_Mead",
                                            optCtrl=list(maxfun = 1e5)))
# GLMER 3 level + slope
form <- TICKS ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
glmer3LevSlope  <- glmer(form, family="poisson",data=grouseticks,
                     control = glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun = 1e5)))

# GLMER 2 level
# data(VerbAgg)
# fmVA <- glmer(r2 ~ Anger + Gender + btype + situ +
#                 (1|id) + (1|item), family = binomial, data =
#                 VerbAgg)

# Sleepstudy
lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Wackier example
data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
lmerSlope2 <- lmer(distance ~ age + (0 + age + nsex|Subject), data=Orthodont)

###############################################
context("Extract Random Effects from merMod")
################################################

test_that("REextract pulls out a data frame", {
  expect_is(REextract(lmerSlope2), "data.frame")
  expect_is(REextract(glmer3LevSlope), "data.frame")
  expect_is(REextract(glmer3Lev), "data.frame")
  expect_is(REextract(lmerSlope1), "data.frame")
})

test_that("REextract issues error with non merMod objects", {
  expect_error(REextract(lm(Reaction ~ Days, sleepstudy)))
  expect_error(REextract(glm(TICKS ~ YEAR + HEIGHT,
                             family="poisson", data=grouseticks)))
})

test_that("REextract gets correct dimensions", {
  expect_equal(ncol(REextract(glmer3Lev)), 3)
  expect_equal(ncol(REextract(lmerSlope1)), 5)
  expect_equal(ncol(REextract(lmerSlope2)), 5)
  expect_equal(ncol(REextract(glmer3LevSlope)), 5)
  expect_equal(nrow(REextract(glmer3Lev)), 584)
  expect_equal(nrow(REextract(lmerSlope1)), 18)
  expect_equal(nrow(REextract(lmerSlope2)), 27)
  expect_equal(nrow(REextract(glmer3LevSlope)), 584)
})

###############################################
context("Fixed effect estimates from posterior")
################################################

test_that("FEsim produces data.frames", {
  expect_is(FEsim(lmerSlope1, nsims=100), "data.frame")
  expect_is(FEsim(lmerSlope2, nsims=100), "data.frame")
  expect_is(FEsim(glmer3Lev, nsims=100), "data.frame")
  expect_is(FEsim(glmer3LevSlope, nsims=100), "data.frame")
})


