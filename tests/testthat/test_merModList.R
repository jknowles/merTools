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
})

test_that("print methods work for merModList", {

})


test_that("ICC function works", {

})
