# test merModList functions

context("Do merModList objects get built and work")

test_that("simple cases work", {
  library(blme)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  split <- sample(x = LETTERS[9:15], size = nrow(d), replace=TRUE)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  out <- split(d, split)
  rm(split)
  g1 <- lmerModList(formula = y~fac1+(1|grp), data=out)
  expect_is(g1, "merModList")
  g2 <- blmerModList(formula = y~fac1+(1|grp), data=out)
  expect_is(g2, "merModList")
  expect_false(class(g1[[1]]) == class(g2[[1]]))

  split <- sample(x = LETTERS[1:20], size = nrow(InstEval), replace=TRUE)
  out <- split(InstEval, split)
  rm(split)
  g1 <- lmerModList(formula = y ~ lectage + studage + (1|d) + (1|dept),
                    data=out)
  expect_is(g1, "merModList")

})

test_that("print methods work for merModList", {
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  split <- sample(x = LETTERS[9:15], size = nrow(d), replace=TRUE)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  out <- split(d, split)
  rm(split)
  g1 <- lmerModList(formula = y~fac1+(1|grp), data=out)
  {
    sink("NUL"); zz <- print(g1); sink()
    }
  expect_null(zz)

})


test_that("ICC function works", {
  ICC1 <- ICC(outcome = "Reaction", group = "Subject", data = sleepstudy)
  expect_is(ICC1, "numeric")
  expect_equal(ICC1, 0.44685, tol = .001)
})
