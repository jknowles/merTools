# # -----------------------------------------------------
# #-------------------------------------------------------
set.seed(101)
library(merTools)

context("Prediction intervals cover for simulated problems")

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
  d$fitted <- predict(g1, d)
  outs <- predictInterval(g1, newdata = d, level = 0.9, n.sims = 500,
                          stat = 'mean', include.resid.var = TRUE)
  outs <- cbind(d, outs); outs$coverage <- FALSE
  outs$coverage <- outs$fitted <= outs$upr & outs$fitted >= outs$lwr
  expect_true(all(outs$coverage))
  expect_less_than(abs(mean(outs$fit - outs$fitted)), .0001)
  expect_less_than(abs(mean(outs$fit - outs$y)), .01)
  rm(outs)
})


test_that("Prediction intervals work for simple GLM example", {
  skip_on_cran()
  set.seed(101)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:50)
  d$y <- simulate(~fac1+(1|grp),family = binomial,
                  newdata=d,
                  newparams=list(beta=c(2,-1,3,-2,1.2), theta=c(.33)))[[1]]
  subD <- d[sample(row.names(d), 1200),]

  g1 <- glmer(y~fac1+(1|grp), data=subD, family = 'binomial')
  d$fitted <- predict(g1, d)
  outs <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 500,
                          stat = 'mean', include.resid.var = FALSE,
                          type = 'linear.prediction')
  outs <- cbind(d, outs); outs$coverage <- FALSE
  outs$coverage <- outs$fitted <= outs$upr & outs$fitted >= outs$lwr
  expect_true(all(outs$coverage))
  expect_less_than(abs(mean(outs$fit - outs$fitted)), .1)
  expect_less_than(abs(mean(outs$fit - outs$y)), 2)

  outs2 <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 500,
                          stat = 'mean', include.resid.var = FALSE,
                          type = 'probability')
  expect_false(identical(outs, outs2))
  expect_true(max(outs2$fit) <= 1)
  expect_true(min(outs2$fit) >= 0)
  expect_true(max(outs2$lwr) <= 1)
  expect_true(min(outs2$lwr) >= 0)
  expect_true(max(outs2$upr) <= 1)
  expect_true(min(outs2$upr) >= 0)
  expect_false(max(outs$fit) <= 1)
  # expect_true(min(outs$fit) < 0)
  expect_false(max(outs$lwr) <= 1)
  expect_false(min(outs$lwr) >= 0)
  expect_false(max(outs$upr) <= 1)
  rm(outs)
})

test_that("Prediction interval respects user input", {
  skip_on_cran()
  set.seed(101)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:25)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  subD <- d[sample(row.names(d), 1000),]

  g1 <- lmer(y~fac1+(1|grp), data=subD)
  d$fitted <- predict(g1, d)
  outs1 <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                          stat = 'mean', include.resid.var = TRUE)
  outs2 <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 500,
                           stat = 'mean', include.resid.var = TRUE)
  outs1a <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 1500,
                           stat = 'mean', include.resid.var = TRUE)
  outs2a <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 1500,
                            stat = 'mean', include.resid.var = TRUE)
  outs3 <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                           stat = 'mean', include.resid.var = FALSE)
  outs3b <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                           stat = 'median', include.resid.var = FALSE)
  expect_more_than(median(outs2$upr - outs1$upr), 0.1)
  expect_more_than(median(outs2a$upr - outs1a$upr), 0.1)
  expect_less_than(median(outs3$upr - outs1$upr), -.2)
  expect_less_than(median(outs3b$upr - outs1a$upr), -.2)
  expect_less_than(mean(outs1$upr - outs1$lwr), mean(outs2$upr - outs2$lwr))
  expect_less_than(mean(outs1$upr - outs1$lwr), mean(outs1a$upr - outs1a$lwr))
  expect_less_than(mean(outs2$upr - outs2$lwr), mean(outs2a$upr - outs2a$lwr))
  expect_false(median(outs3$fit) == median(outs3b$fit))
})


test_that("Predict handles unused and subset of factor levels", {
  skip_on_cran()
  skip_on_travis()
  set.seed(101)
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
  d1 <- InstEval[1:100, ]
  outs1 <- predictInterval(g1, newdata = d1, level = 0.8, n.sims = 500,
                           stat = 'mean', include.resid.var = TRUE)
  d2 <- rbind(d1, InstEval[670:900,])
  outs1a <- predictInterval(g1, newdata = d2, level = 0.8, n.sims = 500,
                            stat = 'mean', include.resid.var=TRUE)[1:100,]
  expect_is(outs1, "data.frame")
  expect_is(outs1a, "data.frame")
  expect_equal(nrow(outs1), 100)
  expect_equal(nrow(outs1a), 100)
  g2 <- lmer(y ~ lectage + studage + (1+lectage|d) + (1|dept), data=InstEval)
  d2 <- InstEval[670:900,]
  outs1a <- predictInterval(g2, newdata = d2, level = 0.8, n.sims = 500,
                            stat = 'mean', include.resid.var=TRUE)
  expect_is(outs1a, "data.frame")
  expect_equal(nrow(outs1a), 231)
})

test_that("Prediction intervals work for multiple parameters per level", {
  skip_on_cran()
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5)))

  outs1 <- predictInterval(glmer3LevSlope, newdata = grouseticks[1:10,])
  expect_is(outs1, "data.frame")
})
