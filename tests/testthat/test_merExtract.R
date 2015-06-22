# -----------------------------------------------------
# Test framework includes tests for multiple intercepts and
# multiple slopes to ensure extraction of random effects
# works in these scenarios
#-------------------------------------------------------

set.seed(51315)
library(lme4)
data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
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
  expect_equal(ncol(REextract(glmer3Lev)), 4)
  expect_equal(ncol(REextract(lmerSlope1)), 6)
  expect_equal(ncol(REextract(lmerSlope2)), 6)
  expect_equal(ncol(REextract(glmer3LevSlope)), 6)
  expect_equal(nrow(REextract(glmer3Lev)), 584)
  expect_equal(nrow(REextract(lmerSlope1)), 18)
  expect_equal(nrow(REextract(lmerSlope2)), 27)
  expect_equal(nrow(REextract(glmer3LevSlope)), 584)
})

# Check names
# Check numerics

###############################################
context("Fixed effect estimates from posterior")
################################################

test_that("FEsim produces data.frames", {
  expect_is(FEsim(lmerSlope1, nsims=100), "data.frame")
  expect_is(FEsim(lmerSlope2, nsims=100), "data.frame")
  expect_is(FEsim(glmer3Lev, nsims=100), "data.frame")
  expect_is(FEsim(glmer3LevSlope, nsims=100), "data.frame")
})

test_that("nsims changes simulation results", {
  expect_false(identical(FEsim(lmerSlope1, nsims = 1000),
                         FEsim(lmerSlope1, nsims = 10)))
})

# numeric checks

###############################################
context("Random effect estimates from posterior")
################################################

test_that("REsim produces data.frames", {
  expect_is(REsim(lmerSlope1, nsims=100), "data.frame")
  expect_is(REsim(lmerSlope2, nsims=100), "data.frame")
  expect_is(REsim(glmer3Lev, nsims=100), "data.frame")
  expect_is(REsim(glmer3LevSlope, nsims=100), "data.frame")
})


###############################################
context("RMSE estimates")
################################################

test_that("RMSE produces correct variable types", {
  expect_is(RMSE.merMod(lmerSlope1), "numeric")
  expect_is(RMSE.merMod(lmerSlope2), "numeric")
  expect_is(RMSE.merMod(lmerSlope1, scale = TRUE), "numeric")
  expect_is(RMSE.merMod(lmerSlope2, scale = TRUE), "numeric")
})


test_that("RMSE respects scale parameter", {
  expect_false(identical(RMSE.merMod(lmerSlope1),
                         RMSE.merMod(lmerSlope1, scale = TRUE)))
  expect_false(identical(RMSE.merMod(lmerSlope2),
                         RMSE.merMod(lmerSlope2, scale = TRUE)))
  expect_less_than(RMSE.merMod(lmerSlope2, scale = TRUE),
                   RMSE.merMod(lmerSlope2))
  expect_less_than(RMSE.merMod(lmerSlope1, scale = TRUE),
                   RMSE.merMod(lmerSlope1))
})
