# Test plotting functions

# Plot functions return gg objects? ----
context("Plot functions return gg objects?")

test_that("Prediction intervals work for simple linear example", {
  skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  subD <- d[sample(row.names(d), 1000),]

  g1 <- lmer(y~fac1+(1|grp), data=subD)
  FE1 <- FEsim(g1)
  p1 <- plotFEsim(FE1)
  expect_is(p1, "gg")
  p1 <- plotREsim(REsim(g1))
  expect_is(p1, "gg")

})
