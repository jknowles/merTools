
set.seed(101)

#Prediction intervals cover for simulated problems----
context("Prediction intervals cover for simulated problems")

test_that("Prediction intervals work for simple linear example", {
  skip_on_travis()
  skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)), seed = 4548)[[1]]
  subD <- d[sample(row.names(d), 1000),]

  g1 <- lmer(y~fac1+(1|grp), data=subD)
  d$fitted <- predict(g1, d)
  #This suppresses the warning about no parallel backend registered
  outs <- suppressWarnings(
    predictInterval(g1, newdata = d, level = 0.9, n.sims = 1000,
                    seed = 468,
                    stat = 'mean', include.resid.var = TRUE)
  )
  outs <- cbind(d, outs); outs$coverage <- FALSE
  outs$coverage <- outs$fitted <= outs$upr & outs$fitted >= outs$lwr
  expect_true(all(outs$coverage))
  expect_lt(abs(mean(outs$fit - outs$fitted)), .0005)
  expect_lt(abs(mean(outs$fit - outs$y)), .01)
  rm(outs)
})


test_that("Prediction intervals work for simple GLMM example", {
  skip_on_travis()
  skip_on_cran()
  set.seed(101)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:50)
  d$y <- simulate(~fac1+(1|grp),family = binomial,
                  newdata=d,
                  newparams=list(beta=c(2,-1,3,-2,1.2), theta=c(.33)),
                  seed =634)[[1]]
  subD <- d[sample(row.names(d), 1200),]

  g1 <- glmer(y~fac1+(1|grp), data=subD, family = 'binomial')
  d$fitted <- predict(g1, d)
  outs <- predictInterval(g1, newdata = d, level = 0.9, n.sims = 500,
                          stat = 'mean', include.resid.var = TRUE,
                          type = 'linear.prediction', seed = 4563)
  outs <- cbind(d, outs); outs$coverage <- FALSE
  outs$coverage <- outs$fitted <= outs$upr & outs$fitted >= outs$lwr
  expect_true(all(outs$coverage))
  expect_lt(abs(mean(outs$fit - outs$fitted)), .1)
  expect_lt(abs(mean(outs$fit - outs$y)), 2)

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
                                 sigma = c(.23)), seed =463)[[1]]
  subD <- d[sample(row.names(d), 1000),]

  g1 <- lmer(y~fac1+(1|grp), data=subD)
  d$fitted <- predict(g1, d)
  outs1 <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                           stat = 'mean', include.resid.var = TRUE, seed=643)
  outs2 <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 500,
                           stat = 'mean', include.resid.var = TRUE, seed=643)
  outs1a <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 1500,
                            stat = 'mean', include.resid.var = TRUE, seed=643)
  outs2a <- predictInterval(g1, newdata = d, level = 0.95, n.sims = 1500,
                            stat = 'mean', include.resid.var = TRUE, seed=643)
  outs3 <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                           stat = 'mean', include.resid.var = FALSE, seed=643)
  outs3b <- predictInterval(g1, newdata = d, level = 0.8, n.sims = 500,
                            stat = 'median', include.resid.var = FALSE, seed=643)
  outs3c <- predictInterval(g1, newdata = d[1, ], level = 0.8, n.sims = 500,
                            stat = 'median', include.resid.var = FALSE, seed=643)

  expect_gt(median(outs2$upr - outs1$upr), 0.1)
  expect_gt(median(outs2a$upr - outs1a$upr), 0.1)
  expect_lt(median(outs3$upr - outs1$upr), -.2)
  expect_lt(median(outs3b$upr - outs1a$upr), -.2)
  expect_lt(mean(outs1$upr - outs1$lwr), mean(outs2$upr - outs2$lwr))
  expect_lt(mean(outs1$upr - outs1$lwr), mean(outs1a$upr - outs1a$lwr))
  expect_lt(mean(outs2$upr - outs2$lwr), mean(outs2a$upr - outs2a$lwr))
  expect_false(median(outs3$fit) == median(outs3b$fit))
  expect_equal(nrow(outs3c), 1)
})

# Prediction works for all combinations of slopes and intercepts----
context("Prediction works for all combinations of slopes and intercepts")

test_that("Predict handles unused and subset of factor levels", {
  skip_on_cran()
  skip_on_travis()
  set.seed(101)
  moddf <- InstEval[sample(rownames(InstEval), 10000), ]
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=moddf)
  d1 <- InstEval[1:100, ]
  outs1 <- suppressWarnings(predictInterval(g1, newdata = d1, level = 0.8, n.sims = 500,
                           stat = 'mean', include.resid.var = TRUE, seed = 4632))
  d2 <- rbind(d1, InstEval[670:900,])
  outs1a <- suppressWarnings(predictInterval(g1, newdata = d2, level = 0.8, n.sims = 500,
                            stat = 'mean', include.resid.var=TRUE, seed = 4632)[1:100,])
  expect_is(outs1, "data.frame")
  expect_is(outs1a, "data.frame")
  expect_equal(nrow(outs1), 100)
  expect_equal(nrow(outs1a), 100)
  g2 <- lmer(y ~ lectage + studage + (1+lectage|d) + (1|dept), data=moddf)
  d2 <- InstEval[670:900,]
  outs1a <- suppressWarnings(predictInterval(g2, newdata = d2, level = 0.8, n.sims = 500,
                            stat = 'mean', include.resid.var=TRUE, seed = 4632))
  expect_is(outs1a, "data.frame")
  expect_equal(nrow(outs1a), 231)
})

rm(list = ls())

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

  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = grouseticks[1:10,]))
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
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
  expect_is(outs1, "data.frame")
  # Message may not be necessary any more
  # expect_message(predictInterval(glmer3LevSlope, newdata = zNew))
})

# Test for new factor levels----
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
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
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
  #Suppressing warning for known degenerate model below
  glmer3LevSlope  <- suppressWarnings(glmer(form, family="binomial",data=grouseticks,
                           control = glmerControl(optimizer="bobyqa",
                                                  optCtrl=list(maxfun = 1e5))))
  zNew <- grouseticks[1:10,]
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:9] <- "100"
  zNew$BROOD[10] <- "101"
  expect_warning(predictInterval(glmer3LevSlope, newdata = zNew),
                 "Currently, predictions for these values are based only on the")
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
  zNew <- grouseticks[1:10,]
  outs2 <- predictInterval(glmer3LevSlope, newdata = zNew)
  expect_is(outs1, "data.frame")
  expect_is(outs2, "data.frame")
  expect_identical(dim(outs1), dim(outs2))
})

# Numeric accuracy----
context("Numeric accuracy")

# Cases
# new factor level for group term

test_that("Median of prediction interval is close to predict.lmer for single group models", {
  set.seed(2311)
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  truPred <- predict(fm1, newdata = sleepstudy)
  newPred <- predictInterval(fm1, newdata = sleepstudy, n.sims = 500,
                             level = 0.9, stat = c("median"),
                             include.resid.var = FALSE, seed = 4563)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/50)

  fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  truPred <- predict(fm1, newdata = sleepstudy)
  newPred <- predictInterval(fm1, newdata = sleepstudy, n.sims = 1500,
                             level = 0.9, stat = c("median"),
                             include.resid.var = FALSE, seed = 9598)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/100)
})

test_that("Median of PI is close to predict.lmer for complex group models", {
  skip_on_cran()
  skip_on_travis()
  set.seed(101)
  moddf <- InstEval[sample(rownames(InstEval), 10000), ]
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=moddf)
  d1 <- moddf[1:200, ]
  newPred <- predictInterval(g1, newdata = d1, level = 0.8, n.sims = 500,
                             stat = 'median', include.resid.var = FALSE,
                             seed = 4563)
  truPred <- predict(g1, newdata = d1)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/100)
  rm(list=ls())
})

test_that("Median of PI is close to predict.glmer for basic and complex grouping", {
  skip_on_cran()
  skip_on_travis()
  set.seed(8496)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:8), fac2 = LETTERS[12:19],
                   obs=1:20)
  d$x <- runif(nrow(d))
  d$y <- simulate(~ x + fac1 + fac2 + (1 + fac1|grp) + (1|obs), family = binomial,
                  newdata=d,
                  newparams=list(beta = rnorm(13),
                                 theta = rnorm(16, 5, 1)), seed = 4563)[[1]]
  subD <- d[sample(row.names(d), 1500),]
# TOO SLOW
  g1 <- glmer(y ~ x + fac1 + fac2 + (1+fac1|grp) + (1|obs), data = subD,
              family = 'binomial',
              control = glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun = 1e5)))
  truPred <- predict(g1, subD)
  newPred <- suppressWarnings(predictInterval(g1, newdata = subD, level = 0.95, n.sims = 2000,
                             stat = 'median', include.resid.var = FALSE,
                             type = 'linear.prediction', seed = 3252))
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/20)
  # This test fails currently
  #   g1 <- glmer(y ~ x +  fac2 + (1 + fac1|grp) + (1|obs), data = subD, family = 'binomial')
  #   truPred <- predict(g1, subD, type = "response")
  #   newPred <- predictInterval(g1, newdata = subD, level = 0.8, n.sims = 500,
  #                              stat = 'median', include.resid.var = FALSE,
  #                              type = 'probability')
  #   expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/20)
  rm(list = ls())
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
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  newPred <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew, level = 0.95,
                             n.sims = 500, stat = 'median',
                             include.resid.var = TRUE, seed = 4563))
  truPred <- predict(glmer3LevSlope, newdata = zNew, allow.new.levels = TRUE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/40)
})

test_that("Prediction intervals work with slope not in fixed effects and data reordered", {
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
  zNew <- grouseticks
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:99] <- "100"
  zNew$BROOD[100] <- "101"
  zNew <- zNew[, c(10, 9, 8, 7, 1, 2, 3, 4, 5, 6, 10)]
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  newPred <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew, level = 0.95,
                             n.sims = 500, stat = 'median',
                             include.resid.var = TRUE, seed = 4563))
  truPred <- predict(glmer3LevSlope, newdata = zNew, allow.new.levels = TRUE)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/20)
})

# Special cases - rank deficiency----
context("Special cases - rank deficiency")

test_that("Prediction intervals are accurate with interaction terms and rank deficiency", {
  skip_on_travis()
  skip_on_cran()
  set.seed(54656)
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

  newPred <- predictInterval(fm, newdata = d2, level = 0.8, n.sims = 1500,
                             stat = 'median', include.resid.var = FALSE,
                             seed = 2342)
  truPred <- predict(fm, newdata = d2)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/15)
  fm2 <- lmer( z ~ a*b + (1+b|r), data=d2)
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  newPred <- suppressWarnings(predictInterval(fm2, newdata = d2, level = 0.8,
                                              n.sims = 1000,
                             stat = 'median', include.resid.var = FALSE))
  truPred <- predict(fm2, newdata = d2)
  expect_is(newPred, "data.frame")
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/10)
})

# Test the simResults----
context("Test the simResults")

test_that("simResults option behaves", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  preds1 <- predictInterval(m1, newdata = sleepstudy[1:5, ])
  preds2 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE)
  expect_null(attr(preds1, "sim.results"))
  expect_is(attr(preds2, "sim.results"), "matrix")
  out <- attr(preds2, "sim.results")
  expect_equal(ncol(out), 1000)
  expect_equal(nrow(out), 5)
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  preds1 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE,
                            which = "random", seed = 23151)
  preds2 <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "fixed",
                            returnSims = TRUE, seed = 23151)
  preds3 <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "all",
                            returnSims = TRUE, seed = 23151)
  preds4 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE, seed = 23151)
  preds1b <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE,
                            which = "random", seed = 23151,
                            include.resid.var = FALSE)
  preds2b <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "fixed",
                            returnSims = TRUE, seed = 23151,
                            include.resid.var = FALSE)
  preds3b <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "all",
                            returnSims = TRUE, seed = 23151,
                            include.resid.var = FALSE)
  preds4b <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE, seed = 23151,
                            include.resid.var = FALSE)
  expect_is(attr(preds1, "sim.results"), "matrix")
  expect_gt(abs(mean(attr(preds1, "sim.results") - attr(preds2, "sim.results"))),
                   200)
  expect_is(attr(preds2, "sim.results"), "matrix")
  expect_gt(abs(mean(attr(preds4, "sim.results") - attr(preds1, "sim.results"))),
                   200)
  expect_gt(abs(mean(attr(preds4, "sim.results") - attr(preds1, "sim.results"))),
                   abs(mean(attr(preds1, "sim.results") - attr(preds2, "sim.results"))))
  expect_is(attr(preds3, "sim.results"), "list")
  expect_is(attr(preds4, "sim.results"), "matrix")
  expect_is(attr(preds1b, "sim.results"), "matrix")
  expect_is(attr(preds2b, "sim.results"), "matrix")
  expect_is(attr(preds3b, "sim.results"), "list")
  expect_is(attr(preds4b, "sim.results"), "matrix")
  # Check that samples are wider for include.resid.var = TRUE
  expect_gt(quantile(attr(preds1, "sim.results"), probs = 0.9) - quantile(attr(preds1b, "sim.results"), probs = 0.9),
                   20)
  expect_lt(quantile(attr(preds1, "sim.results"), probs = 0.1) - quantile(attr(preds1b, "sim.results"), probs = 0.1),
                   -20)
  expect_gt(quantile(attr(preds2, "sim.results"), probs = 0.9) - quantile(attr(preds2b, "sim.results"), probs = 0.9),
                   20)
  expect_lt(quantile(attr(preds2, "sim.results"), probs = 0.1) - quantile(attr(preds2b, "sim.results"), probs = 0.1),
                   -20)
  expect_gt(quantile(attr(preds4, "sim.results"), probs = 0.9) - quantile(attr(preds4b, "sim.results"), probs = 0.9),
                   15)
  expect_lt(quantile(attr(preds4, "sim.results"), probs = 0.1) - quantile(attr(preds4b, "sim.results"), probs = 0.1),
                   -15)
})

# Test out of sample predictions----
context("Test out of sample predictions")

test_that("predictInterval makes predictions without observed outcome", {
  skip_on_travis()
  skip_on_cran()
  possNames <- expand.grid(letters,LETTERS)
  possNames <- paste(possNames[, 1], possNames[, 2])
  newFac <- sample(possNames, 32)
  modData <- data.frame(
    y = rnorm(500),
    x = rnorm(500),
    team_name = sample(newFac, 500, replace = TRUE)
  )
  modData$y[251:500] <- rep(NA, 250)
  m0 <- lmer(y ~ x + (1|team_name), data = modData[1:250,])
  #In the calls below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  testPreds1 <- suppressWarnings(predictInterval(m0, newdata = modData[, c(3, 2, 1)]))
  testPreds2 <- suppressWarnings(predictInterval(m0, newdata = modData[1:250, c(2, 3, 1)]))
  testPreds3 <- suppressWarnings(predictInterval(m0, newdata = modData[251:500,]))
  expect_is(testPreds1, "data.frame")
  expect_is(testPreds2, "data.frame")
  expect_is(testPreds3, "data.frame")
})

# Input validation checks----
context("Input validation checks")


test_that("dplyr objects are successfully coerced", {
  skip_on_cran()
  set.seed(101)
  data(sleepstudy)
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  predData <- sleepstudy %>% group_by(Subject) %>% dplyr::summarise(Days = mean(Days))
  expect_warning(predictInterval(m1, newdata = predData),
                 regexp = "newdata is tbl_df or tbl object from dplyr package", all=FALSE)
  #Suppress the warning that we tested for above
  preds2 <- suppressWarnings(predictInterval(m1, newdata = predData, n.sims=2000))
  expect_is(preds2, "data.frame")
  predData2 <- as.data.frame(predData)
  preds1 <- predictInterval(m1, newdata = predData2, n.sims=2000)
  expect_true(sum(preds1$fit - preds2$fit) > -50 & sum(preds1$fit - preds2$fit) < 50)
})

# Model type warnings for non-binomial GLMM----
context("Model type warnings for non-binomial GLMM")

test_that("Warnings issued", {
  skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:50)
  d$y <- simulate(~fac1+(1|grp),family = poisson,
                  newdata=d,
                  newparams=list(beta=c(2,-1,3,-2,1.2), theta=c(.33)),
                  seed = 5636)[[1]]
  g1 <- glmer(y~fac1+(1|grp), data=d, family = 'poisson')
  expect_warning(predictInterval(g1, newdata = d[1:100,]))
  rm(list = ls())
})

# Test Parallel----
context("Test Parallel")

test_that("parallelization does not throw errors and generates good results", {
  skip_on_cran()
  skip_on_travis()
  library(foreach)
  set.seed(1241)
  #TODO reign in memory usage and cpu time here
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  predA <- predictInterval(m1, newdata = m1@frame, n.sims = 2200, seed = 54,
                           include.resid.var = FALSE, stat = "median")
  predB <- predictInterval(m1, newdata = m1@frame, n.sims = 1750, seed = 54,
                           include.resid.var = FALSE, stat = "median")
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .2)
  predA <- predictInterval(m1, newdata = m1@frame, n.sims = 2500, seed = 2141,
                           include.resid.var = FALSE)
  predB <- predictInterval(m1, newdata = m1@frame, n.sims = 1500, seed = 2141,
                           include.resid.var = FALSE)
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .01)
  moddf <- InstEval[sample(rownames(InstEval), 5000),]
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=moddf)
  predA <- predictInterval(g1, newdata = g1@frame, n.sims = 2500, seed = 2141,
                           include.resid.var = FALSE)
  predB <- predictInterval(g1, newdata = g1@frame, n.sims = 1500, seed = 2141,
                           include.resid.var = FALSE)
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .01)
  predA <- predictInterval(g1, newdata = g1@frame[1:499,], n.sims = 2500, seed = 2141,
                           include.resid.var = TRUE)
  predB <- predictInterval(g1, newdata = g1@frame[1:501,], n.sims = 2500, seed = 2141,
                           include.resid.var = TRUE)
  expect_equal(mean(predA$fit[1:499] - predB$fit[1:499]), 0 , tolerance = .002)
  detach("package:foreach", character.only=TRUE)
})

context("Test returning predict interval components")

# Test the option to return different predictInterval components
m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
data(grouseticks)
grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
grouseticks$YEAR  <- as.numeric(grouseticks$YEAR)
grouseticks$HEIGHT  <- grouseticks$HEIGHT - 462.2
glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks)

test_that("Output is correct dimensions", {
  pred1 <- predictInterval(m1, which = "random")
  pred2 <- predictInterval(m2, which = "fixed")
  pred3 <- predictInterval(m2, which = "all")
  expect_equal(nrow(pred1), nrow(pred2))
  expect_equal(ncol(pred1), 3)
  expect_equal(ncol(pred2), 3)
  expect_equal(ncol(pred3), 5)
  expect_equal(nrow(pred3), 180*3)
  expect_equal(nrow(pred1), 180)
  expect_false(nrow(pred3) == nrow(pred2))
  expect_true(nrow(pred3) > nrow(pred2))
  pred1 <- predictInterval(m1, which = "random", stat = "mean")
  pred2 <- predictInterval(m2, which = "fixed", stat = "mean")
  pred3 <- predictInterval(m2, which = "all", stat = "mean")
  expect_equal(nrow(pred1), nrow(pred2))
  expect_equal(ncol(pred1), 3)
  expect_equal(ncol(pred2), 3)
  expect_equal(ncol(pred3), 5)
  expect_equal(nrow(pred3), 180*3)
  expect_equal(nrow(pred1), 180)
  expect_false(nrow(pred3) == nrow(pred2))
  expect_true(nrow(pred3) > nrow(pred2))
  pred1 <- suppressWarnings(predictInterval(glmer3LevSlope, which = "random",
                                            type = "linear.prediction"))
  pred2 <- suppressWarnings(predictInterval(glmer3LevSlope, which = "random",
                                            type = "probability"))
  pred3 <- suppressWarnings(predictInterval(glmer3LevSlope, which = "all",
                           type = "linear.prediction"))
  pred4 <- suppressWarnings(predictInterval(glmer3LevSlope, which = "all",
                                            type = "probability"))
  expect_equal(nrow(pred1), nrow(pred2))
  expect_equal(nrow(pred3), nrow(pred4))
  expect_equal(ncol(pred1), 3)
  expect_equal(ncol(pred2), 3)
  expect_equal(ncol(pred3), 5)
  expect_equal(ncol(pred3), 5)
  expect_true(mean(pred2$fit) > mean(pred1$fit))
  expect_equal(mean(pred2$fit), 0.5, tol = 0.01)
  expect_equal(mean(pred1$fit), 0.00, tol = 0.05)
  expect_equal(mean(pred4$fit), mean(grouseticks$TICKS_BIN), tol = 0.15)
  expect_false(mean(pred4$fit) == mean(pred3$fit))
  expect_equal(nrow(pred1), 403)
  expect_equal(nrow(pred2), 403)
  expect_equal(nrow(pred3), 403*5)
  expect_equal(nrow(pred4), 403*5)
})



test_that("Compare random, fixed, include-resid", {
  predInt1 <- predictInterval(m1)
  predInt1$width <- predInt1[, 2] - predInt1[, 3]
  predInt2 <- predictInterval(m1, which = "random")
  predInt2$width <- predInt2[, 2] - predInt2[, 3]
  predInt3 <- predictInterval(m1, which = "fixed")
  predInt3$width <- predInt3[, 2] - predInt3[, 3]
  predInt4 <- predictInterval(m1, which = "all")
  predInt4$width <- predInt4[, 3] - predInt4[, 4]
  predInt1b <- predictInterval(m1, include.resid.var = FALSE)
  predInt1b$width <- predInt1b[, 2] - predInt1b[, 3]
  predInt2b <- predictInterval(m1, which = "random", include.resid.var = FALSE)
  predInt2b$width <- predInt2b[, 2] - predInt2b[, 3]
  predInt3b <- predictInterval(m1, which = "fixed", include.resid.var = FALSE)
  predInt3b$width <- predInt3b[, 2] - predInt3b[, 3]
  predInt4b <- predictInterval(m1, which = "all", include.resid.var = FALSE)
  predInt4b$width <- predInt4b[, 3] - predInt4b[, 4]
  # These should be fairly close now since residual variance is the biggest
  expect_true(!all(predInt1$width > predInt2$width))
  expect_true(!all(predInt3$width > predInt2$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt1$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt2$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt3$width))
  #
  expect_true(all(predInt1b$width > predInt2b$width))
  expect_true(all(predInt3b$width != predInt2b$width))
  expect_true(all(predInt4b$width[predInt4b$effect == "combined"] > predInt2b$width))
  expect_true(!all(predInt4b$width[predInt4b$effect == "combined"] > predInt1b$width))
  # Fits
  expect_true(all(predInt1$fit > predInt2$fit))
  expect_true(all(predInt3$fit > predInt2$fit))
  expect_true(all(predInt1$upr > predInt1b$upr))
  expect_true(all(predInt1$lwr < predInt1b$lwr))
  expect_true(all(predInt2$upr > predInt2b$upr))
  expect_true(all(predInt2$lwr < predInt2b$lwr))
  expect_true(all(predInt2$upr > predInt2b$upr))
  expect_true(all(predInt3$lwr < predInt3b$lwr))
  expect_true(all(predInt4$upr > predInt4b$upr))
  expect_true(all(predInt4$lwr < predInt4b$lwr))

  predInt1 <- predictInterval(m2)
  predInt1$width <- predInt1[, 2] - predInt1[, 3]
  predInt2 <- predictInterval(m2, which = "random")
  predInt2$width <- predInt2[, 2] - predInt2[, 3]
  predInt3 <- predictInterval(m2, which = "fixed")
  predInt3$width <- predInt3[, 2] - predInt3[, 3]
  predInt4 <- predictInterval(m2, which = "all")
  predInt4$width <- predInt4[, 3] - predInt4[, 4]
  predInt1b <- predictInterval(m2, include.resid.var = FALSE)
  predInt1b$width <- predInt1b[, 2] - predInt1b[, 3]
  predInt2b <- predictInterval(m2, which = "random", include.resid.var = FALSE)
  predInt2b$width <- predInt2b[, 2] - predInt2b[, 3]
  predInt3b <- predictInterval(m2, which = "fixed", include.resid.var = FALSE)
  predInt3b$width <- predInt3b[, 2] - predInt3b[, 3]
  predInt4b <- predictInterval(m2, which = "all", include.resid.var = FALSE)
  predInt4b$width <- predInt4b[, 3] - predInt4b[, 4]
  # These should be fairly close now since residual variance is the biggest variance
  expect_true(!all(predInt1$width > predInt2$width))
  expect_true(!all(predInt3$width > predInt2$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt1$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt2$width))
  expect_true(!all(predInt4$width[predInt4$effect == "combined"] > predInt3$width))
  #
  expect_true(all(predInt1b$width > predInt2b$width))
  expect_true(all(predInt3b$width != predInt2b$width))
  expect_true(all(predInt4b$width[predInt4b$effect == "combined"] > predInt2b$width))
  expect_true(!all(predInt4b$width[predInt4b$effect == "combined"] > predInt1b$width))
  # Fits
  expect_true(all(predInt1$fit > predInt2$fit))
  expect_true(all(predInt3$fit > predInt2$fit))
  expect_true(all(predInt1$upr > predInt1b$upr))
  expect_true(all(predInt1$lwr < predInt1b$lwr))
  expect_true(all(predInt2$upr > predInt2b$upr))
  expect_true(all(predInt2$lwr < predInt2b$lwr))
  expect_true(all(predInt2$upr > predInt2b$upr))
  expect_true(all(predInt3$lwr < predInt3b$lwr))
  expect_true(all(predInt4$upr > predInt4b$upr))
  expect_true(all(predInt4$lwr < predInt4b$lwr))
})

test_that("Default is set to all effects", {
  predInt1 <- predictInterval(m1, seed = 8231)
  predInt2 <- predictInterval(m1, which = "full", seed = 8231)
  expect_identical(predInt1, predInt2)
  predInt1 <- predictInterval(m1, seed = 8231, include.resid.var = TRUE)
  predInt2 <- predictInterval(m1, which = "full", seed = 8231,
                              include.resid.var = TRUE)
  expect_identical(predInt1, predInt2)
  predInt1 <- predictInterval(m1, seed = 8231, include.resid.var = TRUE)
  predInt2 <- predictInterval(m1, which = "all", seed = 8231,
                              include.resid.var = TRUE)
  expect_equal(mean(predInt1$fit - predInt2$fit[predInt2$effect == "combined"]),
               0, tol = 0.14)
  expect_equal(mean(predInt1$lwr - predInt2$lwr[predInt2$effect == "combined"]),
               0, tol = 0.1)
  expect_equal(mean(predInt1$upr - predInt2$upr[predInt2$effect == "combined"]),
               0, tol=0.1)
})


# Test nested effect specifications----
context("Test nested effect specifications")

test_that("Nested effects can work", {
  skip_on_cran()
  library(ggplot2)
  mod1 <- lmer(sleep_total ~ bodywt + (1|vore/order), data=msleep)
  msleep$combn <- paste(msleep$vore, msleep$order, sep = "__")
  mod2 <- lmer(sleep_total ~ bodywt +  (1|combn) + (1|vore), data=msleep)
  #Suppressing warnings we already tested (coerce tbl and new levels)
  predInt1 <- suppressWarnings(predictInterval(merMod=mod1, newdata=msleep, seed = 548,
                              n.sims = 2000, include.resid.var = FALSE,
                              stat = "median", level = 0.8))
  predInt2 <- suppressWarnings(predictInterval(merMod=mod2, newdata=msleep, seed = 548,
                              n.sims = 2000, include.resid.var = FALSE,
                              stat = "median", level = 0.8))
  expect_is(predInt1, "data.frame")
  expect_is(predInt2, "data.frame")
  expect_equal(mean(predInt1[,1] - predInt2[,1]), 0, tol = sd(predInt1[,1])/20)
  expect_equal(mean(predInt1[,2] - predInt2[,2]), 0, tol = sd(predInt1[,2])/10)
  expect_equal(mean(predInt1[,3] - predInt2[,3]), 0, tol = sd(predInt1[,3])/20)
})


