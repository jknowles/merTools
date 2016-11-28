# Test fastdisplay

set.seed(51315)
library(lme4); library(arm)
# Sleepstudy
lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)


###############################################
context("Fast display")
################################################

test_that("fastdisp pulls out a list", {
  # hack to avoid console output
  {sink("NUL"); zz <- fastdisp(lmerSlope1); sink()}
  expect_is(zz, "list")
  expect_identical(names(zz), c("call", "t.value", "coef", "se",
                                "ngrps", "AIC", "n"))
})

test_that("fastdisp speed is good", {
   {sink("NUL"); t1 <- system.time(force(fastdisp(lmerSlope1)))["elapsed"];
    sink()}
   {sink("NUL"); t2 <- system.time(force(display(lmerSlope1)))["elapsed"];
    sink()}
   expect_lt(t1, t2)
   expect_lt(t1, 0.1)
})
