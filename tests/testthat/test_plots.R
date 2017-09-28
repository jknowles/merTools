# Test plotting functions

library(merTools)
library(testthat)


context("plotFEsim and plotREsim")

test_that("... error functions work properly", {
  require(ggplot2); require(lme4);
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  re1 <- REsim(fm1)
  fe1 <- FEsim(fm1)
  
  # errors for plotREsim
  expect_error(plotREsim(re1, level = -1))
  expect_error(plotREsim(re1, level = 0))
  expect_error(plotREsim(re1, level = 5))
  expect_error(plotREsim(re1, level = "abc"))
  expect_error(plotREsim(re1, level= 0.95, stat= "sd"))
  expect_error(plotREsim(re1, level= 0.95, stat= "abc"))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", sd= -1))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", sd= "TRUE"))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", sigmaScale= -1))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", sigmaScale= "TRUE"))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", oddsRatio= -1))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", oddsRatio= "TRUE"))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", labs= -1))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", labs= "TRUE"))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= NULL))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= list("Subject", "(Intercept)")))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= list(1,2)))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= list(a="Subject", b="(Intercept)")))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= list(groupFctr="Subject", b="(Intercept)")))
  expect_error(plotREsim(re1, level= 0.95, stat= "median", facet= list(a="Subject", term="(Intercept)")))
})


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
  p1 <- plotFEsim(FE1, stat= "median")
  expect_is(p1, "gg")
  p1 <- plotREsim(REsim(g1), stat= "median")
  expect_is(p1, "gg")

})


test_that("... faceting works as designed", {
  skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  subD <- d[sample(row.names(d), 1000),]
  
  g1 <- lmer(y~fac1+(1|grp), data=subD)
  re1 <- REsim(g1)
  p1 <- plotREsim(re1, facet= TRUE)
  p2 <- plotREsim(re1, facet= list(groupFctr= "grp", term= "(Intercept)"))
  
  expect_is(p1, "gg")
  expect_is(p2, "gg")
})