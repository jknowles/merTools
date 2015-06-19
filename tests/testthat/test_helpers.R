# Test helper functions

context("Trimming data frame")



test_that("Trimming results in correct size", {
  data(InstEval)
  trimDat <- trimModelFrame(InstEval)
  expect_more_than(nrow(InstEval), nrow(trimModelFrame(InstEval)))
  expect_equal(nrow(trimDat), 4065)
})


context("subBoot and Theta")

test_that("Can extract theta from a fit model", {
  set.seed(404)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  subD <- d[sample(row.names(d), 1000),]
  g1 <- lmer(y~fac1+(1|grp), data=subD)

  g1b <- lm(y ~ fac1, data = subD)

  expect_equal(thetaExtract(g1), 0.2085, tolerance = .05)
  expect_error(thetaExtract(g1b))

  z1 <- subBoot(g1, 500, FUN = thetaExtract, B = 10)
  expect_is(z1, "data.frame")
  expect_equal(nrow(z1), 11)
  expect_equal(ncol(z1), 2)
})
