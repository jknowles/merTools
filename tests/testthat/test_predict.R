# # -----------------------------------------------------
# #-------------------------------------------------------
set.seed(101)
library(merTools)

context("Prediction intervals cover for simulated problems")

test_that("Prediction intervals work for simple linear example", {
  # skip_on_travis()
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
  skip_on_travis()
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
  skip_on_travis()
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

context("Prediction works for all combinations of slopes and intercepts")

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
  skip_on_travis()
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

test_that("Prediction works for random slopes not in fixed", {
  skip_on_travis()
  skip_on_cran()
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + (1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5)))

  zNew <- grouseticks[1:10,]
  outs1 <- predictInterval(glmer3LevSlope, newdata = zNew)
  expect_is(outs1, "data.frame")
  expect_message(predictInterval(glmer3LevSlope, newdata = zNew))
})


context("Test for new factor levels")

test_that("Prediction intervals work with new factor levels added", {
  skip_on_travis()
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

  zNew <- grouseticks[1:10,]
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:9] <- "100"
  zNew$BROOD[10] <- "101"
  outs1 <- predictInterval(glmer3LevSlope, newdata = zNew)
  expect_is(outs1, "data.frame")
  expect_warning(predictInterval(glmer3LevSlope, newdata = zNew))
})


test_that("Prediction works for factor as a random slope not in fixed", {
  skip_on_travis()
  skip_on_cran()
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ HEIGHT +(1 + YEAR|BROOD) + (1|LOCATION)
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5)))
  zNew <- grouseticks[1:10,]
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:9] <- "100"
  zNew$BROOD[10] <- "101"
  outs1 <- predictInterval(glmer3LevSlope, newdata = zNew)
  zNew <- grouseticks[1:10,]
  outs2 <- predictInterval(glmer3LevSlope, newdata = zNew)
})


context("Numeric accuracy")

# Cases
# new factor level for group term

test_that("Median of prediction interval is close to predict.lmer for single group models", {
  set.seed(2311)
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  truPred <- predict(fm1, newdata = sleepstudy)
  newPred <- predictInterval(fm1, newdata = sleepstudy, n.sims = 500,
                             level = 0.9, stat = c("median"), include.resid.var = FALSE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/50)

  fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  truPred <- predict(fm1, newdata = sleepstudy)
  newPred <- predictInterval(fm1, newdata = sleepstudy, n.sims = 500,
                             level = 0.9, stat = c("median"), include.resid.var = FALSE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/100)
})

test_that("Median of PI is close to predict.lmer for complex group models", {
  skip_on_cran()
  skip_on_travis()
  set.seed(101)
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
  d1 <- InstEval[1:200, ]
  newPred <- predictInterval(g1, newdata = d1, level = 0.8, n.sims = 500,
                           stat = 'median', include.resid.var = FALSE)
  truPred <- predict(g1, newdata = d1)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/100)
})

test_that("Median of PI is close to predict.glmer for basic and complex grouping", {
  skip_on_cran()
  skip_on_travis()
  set.seed(3845)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10), fac2 = LETTERS[10:20],
                   obs=1:25)
  d$x <- runif(nrow(d))
  d$y <- simulate(~ x + fac1 + fac2 + (1 + fac1|grp) + (1|obs), family = binomial,
                  newdata=d,
                  newparams=list(beta = rnorm(16),
                                 theta = rnorm(16, 5, 1)))[[1]]
  subD <- d[sample(row.names(d), 5000),]

  g1 <- glmer(y ~ x + fac1 + fac2 + (1+fac1|grp) + (1|obs), data = subD, family = 'binomial')
  truPred <- predict(g1, subD, type = "response")
  newPred <- predictInterval(g1, newdata = subD, level = 0.95, n.sims = 500,
                          stat = 'median', include.resid.var = FALSE,
                          type = 'probability')
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/40)
  # This test fails currently
#   g1 <- glmer(y ~ x +  fac2 + (1 + fac1|grp) + (1|obs), data = subD, family = 'binomial')
#   truPred <- predict(g1, subD, type = "response")
#   newPred <- predictInterval(g1, newdata = subD, level = 0.8, n.sims = 500,
#                              stat = 'median', include.resid.var = FALSE,
#                              type = 'probability')
#   expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/20)
})

test_that("Prediction intervals work with new factor levels added", {
  skip_on_cran()
  skip_on_travis()
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5)))

  zNew <- grouseticks
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:99] <- "100"
  zNew$BROOD[100] <- "101"
  newPred <- predictInterval(glmer3LevSlope, newdata = zNew, level = 0.95, n.sims = 500,
                           stat = 'median', include.resid.var = TRUE)
  truPred <- predict(glmer3LevSlope, newdata = zNew, allow.new.levels = TRUE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/40)
})

test_that("Prediction intervals work with slope not in fixed effects and data reordered", {
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + (1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5)))
  zNew <- grouseticks
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:99] <- "100"
  zNew$BROOD[100] <- "101"
  zNew <- zNew[, c(10, 9, 8, 7, 1, 2, 3, 4, 5, 6, 10)]
  newPred <- predictInterval(glmer3LevSlope, newdata = zNew, level = 0.95, n.sims = 500,
                             stat = 'median', include.resid.var = TRUE)
  truPred <- predict(glmer3LevSlope, newdata = zNew, allow.new.levels = TRUE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/40)
})

context("Special cases - rank deficiency")

test_that("Prediction intervals work with slope not in fixed effects and data reordered", {
  n <- 20
  x <- y <- rnorm(n)
  z <- rnorm(n)
  r <- sample(1:5, size=n, replace=TRUE)
  d <- data.frame(x,y,z,r)
  d2 <- expand.grid(a=factor(1:4),b=factor(1:4),rep=1:10)
  n <- nrow(d2)
  d2 <- transform(d2,r=sample(1:5, size=n, replace=TRUE),
                  z=rnorm(n))
  d2 <- subset(d2,!(a=="4" & b=="4"))
  fm <- lmer( z ~ a*b + (1|r), data=d2)
  expect_is(predictInterval(fm, newdata = d2[1:10, ]), "data.frame")
})


