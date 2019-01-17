# test subboot

set.seed(51315)
library(lme4)
# Sleepstudy
lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)


###############################################
# Extract theta----
context("Extract theta")
################################################

test_that("extract theta produces a vector", {
  expect_is(thetaExtract(lmerSlope1), "numeric")
  expect_equal(length(thetaExtract(lmerSlope1)), 3)
})

test_that("thetaExtract throws errors for non-merMod objects", {
  expect_error(thetaExtract(lmerSlope1@frame))
  m1 <- lm(mpg ~ disp + hp, data = mtcars)
  expect_error(thetaExtract(m1))
})

###############################################
# subBoot----
context("subBoot")
################################################

test_that("subBoot produces correct output", {
  skip_on_cran()
  out1 <- subBoot(lmerSlope1, n = 100, FUN = thetaExtract, R = 100)
  expect_is(out1, "data.frame")
  expect_equal(ncol(out1), 4)
  expect_equal(nrow(out1), 101)
  out2 <- subBoot(lmerSlope1, n = 100,
                  FUN = function(x) getME(x, "fixef"),
                  R = 100)
  expect_is(out2, "data.frame")
  expect_equal(ncol(out2), 3)
  expect_equal(nrow(out2), 101)
})
