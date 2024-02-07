
#Prediction intervals cover for simulated problems----

test_that("Prediction intervals work for simple linear example", {
  set.seed(51315)
  skip_on_ci()
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
  skip_on_ci()
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
  skip_on_ci()
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
#context("Prediction works for all combinations of slopes and intercepts")

test_that("Predict handles unused and subset of factor levels", {
  skip_on_cran()
  skip_on_ci()
  set.seed(101)
  moddf <- InstEval[sample(rownames(InstEval), 10000), ]
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=moddf)
  d1 <- InstEval[1:100, ]
  outs1 <- suppressWarnings(predictInterval(g1, newdata = d1, level = 0.8, n.sims = 500,
                                            stat = 'mean', include.resid.var = TRUE, seed = 4632))
  d2 <- rbind(d1, InstEval[670:900,])
  outs1a <- suppressWarnings(predictInterval(g1, newdata = d2, level = 0.8, n.sims = 500,
                                             stat = 'mean', include.resid.var=TRUE, seed = 4632)[1:100,])
  expect_s3_class(outs1, "data.frame")
  expect_s3_class(outs1a, "data.frame")
  expect_equal(nrow(outs1), 100)
  expect_equal(nrow(outs1a), 100)

  g2 <- lmer(y ~ lectage + studage + (1+lectage|d) + (1|dept), data=moddf) |>
    suppresWarnings()

  d2 <- InstEval[670:900,]
  outs1a <- suppressWarnings(predictInterval(g2, newdata = d2, level = 0.8, n.sims = 500,
                                             stat = 'mean', include.resid.var=TRUE, seed = 4632))
  expect_s3_class(outs1a, "data.frame")
  expect_equal(nrow(outs1a), 231)
})

rm(list = ls())
