
test_that("Prediction intervals work for multiple parameters per level", {
  skip_on_ci()
  skip_on_cran()
  set.seed(5150)
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  suppressMessages({
    glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                             control = glmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun = 1e5)))
  })

  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = grouseticks[1:10,]))
  expect_s3_class(outs1, "data.frame")
})

test_that("Prediction works for random slopes not in fixed", {
  skip_on_ci()
  skip_on_cran()
  set.seed(5150)
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + (1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  suppressMessages({
    glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                             control = glmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun = 1e5)))
  })

  zNew <- grouseticks[1:10,]
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
  expect_s3_class(outs1, "data.frame")
  # Message may not be necessary any more
  # expect_message(predictInterval(glmer3LevSlope, newdata = zNew))
})

# Test for new factor levels----

test_that("Prediction intervals work with new factor levels added", {
  skip_on_ci()
  skip_on_cran()
  set.seed(5150)
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  suppressMessages({
    glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                             control = glmerControl(optimizer="bobyqa",
                                                    optCtrl=list(maxfun = 1e5)))

  })

  zNew <- grouseticks[1:10,]
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:9] <- "100"
  zNew$BROOD[10] <- "101"
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shouting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
  expect_s3_class(outs1, "data.frame")
  # Test all warnings and messages here
  expect_snapshot(predictInterval(glmer3LevSlope, newdata = zNew))
})


test_that("Prediction works for factor as a random slope not in fixed", {
  skip_on_ci()
  skip_on_cran()
  set.seed(5150)
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  # GLMER 3 level + slope
  form <- TICKS_BIN ~ HEIGHT +(1 + YEAR|BROOD) + (1|LOCATION)
  #Suppressing warning for known degenerate model below
  suppressMessages({
    glmer3LevSlope  <- suppressWarnings(glmer(form, family="binomial",data=grouseticks,
                                              control = glmerControl(optimizer="bobyqa",
                                                                     optCtrl=list(maxfun = 1e5))))
  })
  zNew <- grouseticks[1:10,]
  zNew$BROOD <- as.character(zNew$BROOD)
  zNew$BROOD[1:9] <- "100"
  zNew$BROOD[10] <- "101"

  predictInterval(glmer3LevSlope, newdata = zNew) |>
     expect_warning("Currently, predictions for these values are based only on the") |>
    suppressWarnings()

  outs1 <- suppressWarnings(predictInterval(glmer3LevSlope, newdata = zNew))
  zNew <- grouseticks[1:10,]
  # Expect warnings because we have 0 for our reMeans on location in this model
  # mvtnorm now issues a warning about the indefinite/rank-deficient matrix
  predictInterval(glmer3LevSlope, newdata = zNew) |>
    expect_warning("the matrix is either rank-deficient or not positive definite") |>
    suppressWarnings()
  # Add test to confirm this
  suppressWarnings({
    outs2 <- predictInterval(glmer3LevSlope, newdata = zNew)
  })

  expect_s3_class(outs1, "data.frame")
  expect_s3_class(outs2, "data.frame")
  expect_identical(dim(outs1), dim(outs2))
})

# Numeric accuracy----

# Cases
# new factor level for group term

test_that("Median of prediction interval is close to predict.lmer for single group models", {
  skip_on_cran()
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
  skip_on_ci()
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
  skip_on_ci()
  set.seed(8496)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:8), fac2 = LETTERS[12:19],
                   obs=1:20)
  d$x <- runif(nrow(d))
  suppressMessages({
    d$y <- simulate(~ x + fac1 + fac2 + (1 + fac1|grp) + (1|obs), family = binomial,
                    newdata=d,
                    newparams=list(beta = rnorm(13),
                                   theta = rnorm(16, 5, 1)), seed = 4563)[[1]]
  })
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
  skip_on_ci()
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
  skip_on_ci()
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
