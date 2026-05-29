# Regression tests for #124:
# https://github.com/jknowles/merTools/issues/124
#
# For a nested model `y ~ 1 + (1|a/b)`, a prediction frame that mixes observed
# and unobserved levels of the *interaction* grouping factor (`a:b`) used to
# produce non-reproducible results. The mapping from the model matrix back to a
# grouping level used `max.col()`, which (a) returns a column even for an
# all-zero row (an unobserved interaction level) and (b) breaks the resulting
# tie *at random*. So an unobserved `a:b` level silently borrowed a randomly
# chosen observed level's random effect, making the point estimate depend on the
# RNG seed and differ between batch vs. row-by-row prediction.
#
# The documented, correct behavior is that an unobserved grouping level
# contributes zero from that random-effect term (predictions fall back to the
# fixed effects plus any *observed* higher-level random effects, plus residual
# variation). These tests pin that behavior.

make_nested_fit <- function() {
  training_df <- data.frame(
    class = c("c1", "c1", "c1", "c1", "c1", "c1", "c1", "c1", "c2"),
    pupil = c("p1", "p1", "p1", "p2", "p3", "p3", "p4", "p5", "p7"),
    grade = c(4.5, 4.8, 4.7, 5.0, 4.1, 4.0, 3.5, 5.5, 5.5)
  )
  suppressMessages(suppressWarnings(
    lme4::lmer(grade ~ 1 + (1 | class / pupil), data = training_df)
  ))
}

test_that("predictInterval is reproducible for partially-observed nested levels (#124)", {
  model <- make_nested_fit()
  # Row 3 (class = c1, pupil = p6) is unobserved at the class:pupil level.
  query_df <- data.frame(
    class = c("c1", "c1", "c1"),
    pupil = c("p1", "p2", "p6")
  )

  s1 <- suppressWarnings(
    predictInterval(model, query_df, seed = 1, n.sims = 4000)
  )
  s2 <- suppressWarnings(
    predictInterval(model, query_df, seed = 2, n.sims = 4000)
  )

  # The unobserved-level row's point estimate must be stable across seeds.
  # Before the fix this swung by ~0.3 (random column selection); after the fix
  # it varies only by ordinary Monte Carlo noise (< 0.1 here).
  expect_lt(abs(s1$fit[3] - s2$fit[3]), 0.1)

  # The genuinely observed rows were always stable; confirm they still are.
  expect_lt(abs(s1$fit[1] - s2$fit[1]), 0.1)
  expect_lt(abs(s1$fit[2] - s2$fit[2]), 0.1)
})

test_that("unobserved nested level falls back to fixed + observed higher-level RE (#124)", {
  model <- make_nested_fit()
  query_df <- data.frame(
    class = c("c1", "c1", "c1"),
    pupil = c("p1", "p2", "p6")
  )

  pred <- suppressWarnings(
    predictInterval(model, query_df, seed = 1, n.sims = 4000)
  )

  # For an unobserved class:pupil level, the class:pupil contribution is zero,
  # so the fitted value should equal the intercept plus the (observed) class
  # random effect for c1 -- it must NOT borrow a random pupil effect.
  expected_fit <- unname(
    lme4::fixef(model)[["(Intercept)"]] +
      lme4::ranef(model)$class["c1", "(Intercept)"]
  )
  expect_equal(pred$fit[3], expected_fit, tolerance = 0.1)
})

test_that("nested-level prediction is consistent batch vs. row-by-row (#124)", {
  model <- make_nested_fit()
  query_df <- data.frame(
    class = c("c1", "c1", "c1"),
    pupil = c("p1", "p2", "p6")
  )

  batch <- suppressWarnings(
    predictInterval(model, query_df, seed = 1, n.sims = 4000)
  )
  alone <- suppressWarnings(
    predictInterval(model, query_df[3, , drop = FALSE], seed = 1, n.sims = 4000)
  )

  # The unobserved-level row should yield essentially the same estimate whether
  # predicted in a batch with observed rows or on its own.
  expect_lt(abs(batch$fit[3] - alone$fit[1]), 0.1)
})
