set.seed(11213)

test_that("fastdisp.merMod returns a list invisibly with expected components", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  result <- fastdisp(m1)
  expect_type(result, "list")
  expect_true(all(c("call", "coef", "se", "ngrps", "AIC", "n") %in% names(result)))
  expect_type(result$coef, "double")
  expect_type(result$se, "double")
  expect_true(is.numeric(result$AIC))
  expect_true(is.numeric(result$n))
})

test_that("fastdisp.merMod return value can be captured and is non-null", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  result <- fastdisp(m1)
  expect_false(is.null(result))
  expect_equal(names(result$coef), c("(Intercept)", "Days"))
})

test_that("fastdisp.merModList returns a list with expected components", {
  skip_on_cran()
  sim_list <- replicate(
    n = 3,
    expr = sleepstudy[sample(row.names(sleepstudy), 180), ],
    simplify = FALSE
  )
  ml <- lmerModList("Reaction ~ Days + (1 | Subject)", data = sim_list)
  result <- fastdisp(ml)
  expect_type(result, "list")
  expect_true(all(c("call", "coef", "se", "ngrps", "AIC", "n") %in% names(result)))
  expect_false(is.null(result))
})
