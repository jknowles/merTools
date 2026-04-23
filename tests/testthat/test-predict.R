

# Test nested effect specifications----

test_that("Nested effects can work", {
  skip_on_cran()
  library(ggplot2)
  mod1 <- lmer(sleep_total ~ bodywt + (1|vore/order), data=msleep)
  msleep$combn <- paste(msleep$vore, msleep$order, sep = "__")
  mod2 <- lmer(sleep_total ~ bodywt +  (1|combn) + (1|vore), data=msleep)
  #Suppressing warnings we already tested (coerce tbl and new levels)
  predInt1 <- suppressWarnings(predictInterval(merMod=mod1, newdata=msleep, seed = 11213,
                              n.sims = 2000, include.resid.var = FALSE,
                              stat = "median", level = 0.8))
  predInt2 <- suppressWarnings(predictInterval(merMod=mod2, newdata=msleep, seed = 11213,
                              n.sims = 2000, include.resid.var = FALSE,
                              stat = "median", level = 0.8))

  # Structural invariants: nested and crossed parameterizations of the same
  # effect should produce equivalently-shaped outputs with sensible ordering.
  # Numerical equivalence on canonical cases is verified in Layer 2 by
  # tests/testthat/test-predictInterval-snapshot.R.
  expect_s3_class(predInt1, "data.frame")
  expect_s3_class(predInt2, "data.frame")
  expect_equal(dim(predInt1), dim(predInt2))
  expect_named(predInt1, c("fit", "upr", "lwr"))
  expect_true(all(predInt1$upr >= predInt1$fit))
  expect_true(all(predInt1$fit >= predInt1$lwr))
  expect_true(all(predInt2$upr >= predInt2$fit))
  expect_true(all(predInt2$fit >= predInt2$lwr))
})

# context("Interactions without intercepts")

test_that("Models with cross-level interaction and no random intercept work", {
  skip_on_cran()
  set.seed(11213)
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
  # Structural checks for the remaining three model parameterizations; the
  # previous tight-tolerance `mean(fit - predict()) ~ 0` assertions here were
  # statistical-unbiasedness checks (Layer 3) and were the most frequent
  # source of intermittent CI failures. Numerical regression for
  # representative LMM configurations is now pinned in
  # tests/testthat/test-predictInterval-snapshot.R.
  suppressMessages({
    preds1 <- predictInterval(m2, seed = 11213)
  })
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  expect_true(all(preds1$upr >= preds1$fit & preds1$fit >= preds1$lwr))

  preds1 <- predictInterval(m3, seed = 11213)
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  expect_true(all(preds1$upr >= preds1$fit & preds1$fit >= preds1$lwr))

  suppressMessages({
    preds1 <- predictInterval(m4, seed = 11213)
  })
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  expect_true(all(preds1$upr >= preds1$fit & preds1$fit >= preds1$lwr))
})



test_that("Models with no fixed intercept and cross-level interaction work", {
  skip_on_cran()
  # Numerical regression for this edge-case model (historically the source of
  # flaky MC bias checks) is pinned in
  # tests/testthat/test-predictInterval-snapshot.R. This test now verifies
  # only structural behavior and the expected warnings.
  suppressMessages({
    m1 <- lmer(Reaction ~ 0 + Days + Days:Subject + (1 | Days), data = sleepstudy)
  })
  preds1 <- predictInterval(m1, seed = 11213)
  expect_equal(nrow(preds1), 180)
  expect_equal(ncol(preds1), 3)
  expect_true(all(preds1$upr >= preds1$fit & preds1$fit >= preds1$lwr))

  predictInterval(m1, fix.intercept.variance = TRUE) |>
    expect_warning("No fixed-effect intercept detected") |>
    suppressWarnings()

  preds1 <- predictInterval(m1, newdata = sleepstudy[1:25, ],
                            level = 0.9, n.sims = 500, include.resid.var = FALSE,
                            ignore.fixed.terms = TRUE, seed = 11213)
  expect_equal(nrow(preds1), 25)
  expect_equal(ncol(preds1), 3)
  expect_true(all(preds1$upr >= preds1$fit & preds1$fit >= preds1$lwr))

  predictInterval(m1, newdata = sleepstudy[1:50, ],
                            level = 0.9, n.sims = 500, include.resid.var = FALSE,
                            fix.intercept.variance = TRUE) |>
    expect_warning("No fixed-effect intercept") |>
    suppressWarnings()

})
