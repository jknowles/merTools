# Tests for the new.levels argument of predictInterval().
#
# For a grouping level not present in the fitted model, the historical (and
# default) behavior is new.levels = "zero": the random effect is dropped and the
# prediction rests on the fixed effects plus residual variation. new.levels =
# "draw" instead samples that group's effect from the estimated random-effect
# covariance (VarCorr), so the interval reflects between-group uncertainty --
# the analogue of brms::posterior_predict(allow_new_levels = TRUE).

make_slope_fit <- function() {
  suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  ))
}

test_that("new.levels defaults to 'zero' (backward compatible)", {
  m <- make_slope_fit()
  nd <- data.frame(Days = 0:4, Subject = "NEW_SUBJECT")
  def  <- suppressWarnings(predictInterval(m, nd, level = 0.9, n.sims = 2000, seed = 1))
  zero <- suppressWarnings(predictInterval(m, nd, level = 0.9, n.sims = 2000, seed = 1,
                                           new.levels = "zero"))
  expect_identical(def, zero)
})

test_that("new.levels='draw' is wider than 'zero' for an unobserved group", {
  m <- make_slope_fit()
  nd <- data.frame(Days = 0:4, Subject = "NEW_SUBJECT")
  zero <- suppressWarnings(predictInterval(m, nd, level = 0.9, n.sims = 4000, seed = 1,
                                           new.levels = "zero"))
  draw <- suppressWarnings(predictInterval(m, nd, level = 0.9, n.sims = 4000, seed = 1,
                                           new.levels = "draw"))
  expect_gt(mean(draw$upr - draw$lwr), mean(zero$upr - zero$lwr))
})

test_that("new.levels='draw' width matches the theoretical predictive SD", {
  m <- make_slope_fit()
  nd <- data.frame(Days = 0:4, Subject = "NEW_SUBJECT")
  draw <- suppressWarnings(predictInterval(m, nd, level = 0.9, n.sims = 8000, seed = 1,
                                           new.levels = "draw"))
  Sigma <- as.matrix(lme4::VarCorr(m)$Subject)        # RE covariance (intercept, Days)
  sig2  <- sigma(m)^2
  z <- qnorm(0.95)
  # for each Days value, total predictive variance = z' Sigma z + residual
  expected_width <- vapply(nd$Days, function(d) {
    zvec <- c(1, d)
    2 * z * sqrt(as.numeric(t(zvec) %*% Sigma %*% zvec) + sig2)
  }, numeric(1))
  observed_width <- draw$upr - draw$lwr
  # Monte Carlo, so allow a modest relative tolerance
  expect_equal(observed_width, expected_width, tolerance = 0.08)
})

test_that("observations sharing a new level share the sampled effect", {
  m <- make_slope_fit()
  # two rows, same new subject, same Days -> identical sampled random effect.
  # Exclude residual variance so only the (shared) RE component remains; with
  # residual noise on, each row gets an independent residual draw and they would
  # differ even though the underlying RE draw is shared.
  nd <- data.frame(Days = c(3, 3), Subject = c("NEW_A", "NEW_A"))
  pi <- suppressWarnings(predictInterval(m, nd, which = "random", level = 0.9,
                                         n.sims = 1000, seed = 1, new.levels = "draw",
                                         include.resid.var = FALSE, returnSims = TRUE))
  sims <- attr(pi, "sim.results")
  if (nrow(sims) == nrow(nd)) sims <- t(sims)
  # same new level + same covariates -> the two rows' draws coincide
  expect_equal(sims[, 1], sims[, 2])
})
