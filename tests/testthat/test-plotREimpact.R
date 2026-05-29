# Tests for plotREimpact() -- overlaying one or more REimpact() results (#84)

test_that("plotREimpact returns a ggplot for a single REimpact result", {
  skip_on_cran()
  set.seed(11213)
  m1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  imp <- suppressWarnings(
    REimpact(m1, newdata = sleepstudy[1, ], groupFctr = "Subject",
             breaks = 4, n.sims = 100)
  )
  p <- plotREimpact(imp)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plotREimpact overlays a named list of REimpact results (#84)", {
  skip_on_cran()
  set.seed(11213)
  m1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  a <- suppressWarnings(
    REimpact(m1, newdata = sleepstudy[1, ], groupFctr = "Subject",
             breaks = 4, n.sims = 100)
  )
  b <- suppressWarnings(
    REimpact(m1, newdata = sleepstudy[2, ], groupFctr = "Subject",
             breaks = 4, n.sims = 100)
  )
  p <- plotREimpact(list("Case A" = a, "Case B" = b))
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  # The two series must both be represented in the built plot data.
  expect_true(any(grepl("colour|color", names(built$data[[length(built$data)]]))))
})

test_that("plotREimpact honors facet and level arguments", {
  skip_on_cran()
  set.seed(11213)
  m1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  imp <- suppressWarnings(
    REimpact(m1, newdata = sleepstudy[1:2, ], groupFctr = "Subject",
             breaks = 3, n.sims = 100)
  )
  p_wide <- plotREimpact(imp, level = 0.99)
  p_narrow <- plotREimpact(imp, level = 0.80)
  d_wide <- ggplot2::ggplot_build(p_wide)$data[[1]]
  d_narrow <- ggplot2::ggplot_build(p_narrow)$data[[1]]
  # A wider confidence level produces wider-or-equal intervals.
  expect_true(all((d_wide$ymax - d_wide$ymin) >= (d_narrow$ymax - d_narrow$ymin) - 1e-8))
  expect_s3_class(plotREimpact(imp, facet = FALSE), "ggplot")
})

test_that("plotREimpact validates its inputs", {
  expect_error(plotREimpact(data.frame(a = 1)), "REimpact")
  expect_error(plotREimpact(42), "data.frame")
  set.seed(11213)
  m1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  imp <- suppressWarnings(
    REimpact(m1, newdata = sleepstudy[1, ], groupFctr = "Subject",
             breaks = 3, n.sims = 100)
  )
  expect_error(plotREimpact(imp, level = 1.5), "between 0 and 1")
  expect_error(plotREimpact(imp, level = 0), "between 0 and 1")
})
