# test subboot

set.seed(51315)
library(lme4)
# Sleepstudy
lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)


###############################################
# Extract theta----
context("Extract theta using subBoot")
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
  # Subbooot returns errors here
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

context("subBoot glmer models")



test_that("subBoot produces correct glmer output", {
  skip_on_cran()
  d <- expand.grid(fac1 = LETTERS[1:5],
                   grp = letters[11:20],
                   obs = 1:50)
  d$y <- simulate(~fac1 + (1 | grp), family = binomial,
                  newdata = d,
                  newparams = list( beta = c(2,-1,3,-2,1.2),
                                    theta = c(.33)),
                  seed =634)[[1]]
  subD <- d[sample(row.names(d), 1200), ]

  g1 <- glmer(y~fac1+(1|grp), data=subD, family = 'binomial')

  out1 <- subBoot(g1, n = 1000, FUN = thetaExtract, R = 10)
  expect_is(out1, "data.frame")
  expect_equal(ncol(out1), 2)
  expect_equal(nrow(out1), 11)
  #
  out2 <- subBoot(g1, n = 500,
                  FUN = function(x) getME(x, "fixef"),
                  R = 10)
  expect_is(out2, "data.frame")
  expect_equal(ncol(out2), 6)
  expect_equal(nrow(out2), 11)

})
