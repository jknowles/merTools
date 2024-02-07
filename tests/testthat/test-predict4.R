
# Test Parallel----


test_that("parallelization does not throw errors and generates good results", {
  skip_on_cran()
  skip_on_ci()
  library(foreach)
  set.seed(1241)
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  predA <- predictInterval(m1, newdata = m1@frame, n.sims = 2200, seed = 54,
                           include.resid.var = FALSE, stat = "median") |> suppressWarnings()
  predB <- predictInterval(m1, newdata = m1@frame, n.sims = 1750, seed = 54,
                           include.resid.var = FALSE, stat = "median")
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .2)
  predA <- predictInterval(m1, newdata = m1@frame, n.sims = 2500, seed = 2141,
                           include.resid.var = FALSE)
  predB <- predictInterval(m1, newdata = m1@frame, n.sims = 1500, seed = 2141,
                           include.resid.var = FALSE)
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .01)
  moddf <- InstEval[sample(rownames(InstEval), 5000),]
  g1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data = moddf)
  predA <- predictInterval(g1, newdata = g1@frame, n.sims = 2500, seed = 2141,
                           include.resid.var = FALSE)
  predB <- predictInterval(g1, newdata = g1@frame, n.sims = 1500, seed = 2141,
                           include.resid.var = FALSE)
  expect_equal(mean(predA$fit - predB$fit), 0 , tolerance = .01)
  predA <- predictInterval(g1, newdata = g1@frame[1:499,], n.sims = 2500, seed = 2141,
                           include.resid.var = TRUE)
  predB <- predictInterval(g1, newdata = g1@frame[1:501,], n.sims = 2500, seed = 2141,
                           include.resid.var = TRUE)
  expect_equal(mean(predA$fit[1:499] - predB$fit[1:499]), 0 , tolerance = .0025)
  detach("package:foreach", character.only=TRUE)
})

#context("Test returning predict interval components")


test_that("Output is correct dimensions", {
  skip_on_cran()
  ###########################################
  # Test the option to return different predictInterval components
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  data(grouseticks)
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  grouseticks$YEAR  <- as.numeric(grouseticks$YEAR)
  grouseticks$HEIGHT  <- grouseticks$HEIGHT - 462.2
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks) |>
    suppressMessages()

  #############################################
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
  expect_equal(mean(pred2$fit), 0.5, tolerance = 0.05)
  expect_equal(mean(pred1$fit), 0.00, tolerance = 0.05)
  # Tolerance meaning has changed
  expect_equal(mean(pred4$fit), mean(grouseticks$TICKS_BIN), tolerance = 0.2)
  expect_false(mean(pred4$fit) == mean(pred3$fit))
  expect_equal(nrow(pred1), 403)
  expect_equal(nrow(pred2), 403)
  expect_equal(nrow(pred3), 403*5)
  expect_equal(nrow(pred4), 403*5)
})



test_that("Compare random, fixed, include-resid", {
  skip_on_cran()
  ######################################
  # Test the option to return different predictInterval components
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  data(grouseticks)
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  grouseticks$YEAR  <- as.numeric(grouseticks$YEAR)
  grouseticks$HEIGHT  <- grouseticks$HEIGHT - 462.2
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks) |>
    suppressMessages()
  ################################################
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
  skip_on_cran()
  ######################################
  # Test the option to return different predictInterval components
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  form <- TICKS_BIN ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
  data(grouseticks)
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  grouseticks$YEAR  <- as.numeric(grouseticks$YEAR)
  grouseticks$HEIGHT  <- grouseticks$HEIGHT - 462.2
  glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks) |>
    suppressMessages()
  ################################################
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
               0, tolerance =0.14)
  expect_equal(mean(predInt1$lwr - predInt2$lwr[predInt2$effect == "combined"]),
               0, tolerance =0.1)
  expect_equal(mean(predInt1$upr - predInt2$upr[predInt2$effect == "combined"]),
               0, tolerance=0.1)
})
