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
# Extract Random Effects from merMod----
################################################

test_that("REextract pulls out a data frame", {
  expect_s3_class(REextract(lmerSlope2), "data.frame")
  expect_s3_class(REextract(glmer3LevSlope), "data.frame")
  expect_s3_class(REextract(glmer3Lev), "data.frame")
  expect_s3_class(REextract(lmerSlope1), "data.frame")
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
# Fixed effect estimates from posterior----

################################################

test_that("FEsim produces data.frames", {
  expect_s3_class(FEsim(lmerSlope1, n.sims=100), "data.frame")
  expect_s3_class(FEsim(lmerSlope2, n.sims=100), "data.frame")
  expect_s3_class(FEsim(glmer3Lev, n.sims=100), "data.frame")
  expect_s3_class(FEsim(glmer3LevSlope, n.sims=100), "data.frame")
})

test_that("n.sims changes simulation results", {
  expect_false(identical(FEsim(lmerSlope1, n.sims = 1000),
                         FEsim(lmerSlope1, n.sims = 10)))
})

# numeric checks

###############################################
# Random effect estimates from posterior----
################################################

test_that("REsim produces data.frames", {
  expect_s3_class(REsim(lmerSlope1, n.sims=100), "data.frame")
  expect_s3_class(REsim(lmerSlope2, n.sims=100), "data.frame")
  expect_s3_class(REsim(glmer3Lev, n.sims=100), "data.frame")
  expect_s3_class(REsim(glmer3LevSlope, n.sims=100), "data.frame")
})


###############################################
# RMSE estimates----
################################################

test_that("RMSE produces correct variable types", {
  expect_type(RMSE.merMod(lmerSlope1), "double")
  expect_type(RMSE.merMod(lmerSlope2), "double")
  expect_type(RMSE.merMod(lmerSlope1, scale = TRUE), "double")
  expect_type(RMSE.merMod(lmerSlope2, scale = TRUE), "double")
})


test_that("RMSE respects scale parameter", {
  expect_false(identical(RMSE.merMod(lmerSlope1),
                         RMSE.merMod(lmerSlope1, scale = TRUE)))
  expect_false(identical(RMSE.merMod(lmerSlope2),
                         RMSE.merMod(lmerSlope2, scale = TRUE)))
  expect_lt(RMSE.merMod(lmerSlope2, scale = TRUE),
                   RMSE.merMod(lmerSlope2))
  expect_lt(RMSE.merMod(lmerSlope1, scale = TRUE),
                   RMSE.merMod(lmerSlope1))
})
