

# Test nested effect specifications----

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
  expect_s3_class(predInt1, "data.frame")
  expect_s3_class(predInt2, "data.frame")
  expect_equal(mean(predInt1[,1] - predInt2[,1]), 0, tolerance = sd(predInt1[,1])/20)
  expect_equal(mean(predInt1[,2] - predInt2[,2]), 0, tolerance = sd(predInt1[,2])/10)
  expect_equal(mean(predInt1[,3] - predInt2[,3]), 0, tolerance = sd(predInt1[,3])/20)
})

# context("Interactions without intercepts")

test_that("Models with cross-level interaction and no random intercept work", {
  skip_on_cran()
  #################################
  sleepstudy$Test <- rep(sample(c(TRUE, FALSE), length(unique(sleepstudy$Subject)),
                                replace = TRUE), each = 10)
  m1 <- lmer(Reaction ~ Days:Test + (0 + Days | Subject), data = sleepstudy)

  sleepstudy$cars <- sleepstudy$Days*3
  m2 <- lmer(Reaction ~ cars:Test + (0 + Days | Subject), data = sleepstudy)
  m3 <- lmer(Reaction ~ cars:Test + (1 | Subject), data = sleepstudy)
  m4 <- lmer(Reaction ~ cars:Test + (0 + cars | Subject), data = sleepstudy)

  ###################################
  suppressMessages({
    preds1 <- predictInterval(m1)
  })
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  expect_message(predictInterval(m1))
  suppressMessages({
    preds1 <- predictInterval(m1, newdata = sleepstudy[1:10, ],
                              level = 0.9, n.sims = 500, include.resid.var = FALSE,
                              fix.intercept.variance = TRUE)
  })
  expect_equal(nrow(preds1), 10)
  expect_equal(ncol(preds1), 3)
  suppressMessages({
    preds1 <- predictInterval(m1, newdata = sleepstudy[1:10, ],
                              level = 0.9, n.sims = 500, include.resid.var = FALSE,
                              ignore.fixed.terms = TRUE)
  })

  expect_equal(nrow(preds1), 10)
  expect_equal(ncol(preds1), 3)
  suppressMessages({
    preds2 <- predictInterval(m1, newdata = sleepstudy[1:10, ],
                              level = 0.9, n.sims = 500, include.resid.var = FALSE,
                              ignore.fixed.terms = FALSE)
  })
  expect_equal(nrow(preds2), 10)
  expect_equal(ncol(preds2), 3)
  expect_false(any(preds1$fit == preds2$fit))
  rm(preds1, preds2)
  suppressMessages({
    preds1 <- predictInterval(m2)
  })

  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  truPred <- predict(m2)
  expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/100)

  #
  preds1 <- predictInterval(m3)
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  truPred <- predict(m3)
  expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/100)
  #
  suppressMessages({
    preds1 <- predictInterval(m4)
  })

  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  truPred <- predict(m4)
  expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/100)
})



test_that("Models with cross-level interaction and no random intercept work", {
  skip_on_cran()
  suppressMessages({
    m1 <- lmer(Reaction ~ 0 + Days + Days:Subject + (1 | Days), data = sleepstudy)
  })
  preds1 <- predictInterval(m1)
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  predictInterval(m1, fix.intercept.variance = TRUE) |>
    expect_warning("No fixed-effect intercept detected") |>
    suppressWarnings()

  preds1 <- predictInterval(m1, newdata = sleepstudy[1:10, ],
                            level = 0.9, n.sims = 500, include.resid.var = FALSE,
                            ignore.fixed.terms = TRUE)
  expect_equal(nrow(preds1), 10)
  expect_equal(ncol(preds1), 3)
  truPred <- predict(m1, newdata = sleepstudy[1:10,])
  expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/100)

  # This is less close
  preds1 <- predictInterval(m1, newdata = sleepstudy[1:50, ],
                            level = 0.9, n.sims = 500, include.resid.var = FALSE,
                            ignore.fixed.terms = FALSE)
  expect_equal(nrow(preds1), 50)
  expect_equal(ncol(preds1), 3)
  truPred <- predict(m1, newdata = sleepstudy[1:50,])
  expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/25)


  predictInterval(m1, newdata = sleepstudy[1:50, ],
                            level = 0.9, n.sims = 500, include.resid.var = FALSE,
                            fix.intercept.variance = TRUE) |>
    expect_warning("No fixed-effect intercept") |>
    suppressWarnings()

    expect_failure({
    expect_equal(mean(preds1$fit - truPred), 0, tolerance = sd(truPred)/100)
  })


})
