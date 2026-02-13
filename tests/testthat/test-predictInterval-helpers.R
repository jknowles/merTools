# Tests for predictInterval helper functions

# Test simulate_residual_variance() ----

test_that("simulate_residual_variance returns correct length for LMM", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  set.seed(123)
  sigma_vec <- merTools:::simulate_residual_variance(m1, n.sims = 100)
  expect_length(sigma_vec, 100)
  expect_true(all(sigma_vec > 0))
})

test_that("simulate_residual_variance respects seed set externally", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  set.seed(123)
  sigma1 <- merTools:::simulate_residual_variance(m1, n.sims = 100)
  set.seed(123)
  sigma2 <- merTools:::simulate_residual_variance(m1, n.sims = 100)
  set.seed(456)
  sigma3 <- merTools:::simulate_residual_variance(m1, n.sims = 100)
  expect_identical(sigma1, sigma2)
  expect_false(identical(sigma1, sigma3))
})

test_that("simulate_residual_variance returns NULL for binomial GLMM", {
   skip_on_cran()
   d1 <- cbpp
   d1$y <- d1$incidence / d1$size
   gm1 <- glmer(y ~ period + (1 | herd), family = binomial, data = d1, weights = d1$size)
   sigma_vec <- merTools:::simulate_residual_variance(gm1, n.sims = 100)
   expect_null(sigma_vec)
 })

test_that("simulate_residual_variance warns for non-binomial GLMM", {
  skip_on_cran()
  skip_on_ci()
  d <- data.frame(
    y = rpois(100, 3),
    x = rnorm(100),
    g = factor(rep(1:10, each = 10))
  )
  gm1 <- glmer(y ~ x + (1 | g), family = poisson, data = d)
  expect_warning(
    merTools:::simulate_residual_variance(gm1, n.sims = 100),
    "not tested"
  )
})


# Test simulate_fixed_effects() ----

test_that("simulate_fixed_effects returns correct dimensions", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  fixed_mat <- merTools:::simulate_fixed_effects(
    m1,
    newdata = sleepstudy[1:10, ],
    n.sims = 100
  )
  expect_equal(nrow(fixed_mat), 10)
  expect_equal(ncol(fixed_mat), 100)
})

test_that("simulate_fixed_effects respects ignore.fixed.terms", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  set.seed(123)
  fixed1 <- merTools:::simulate_fixed_effects(
    m1,
    newdata = sleepstudy[1:5, ],
    n.sims = 500,
    ignore.fixed.terms = NULL
  )
  set.seed(123)
  fixed2 <- merTools:::simulate_fixed_effects(
    m1,
    newdata = sleepstudy[1:5, ],
    n.sims = 500,
    ignore.fixed.terms = "(Intercept)"
  )
  # With ignore.fixed.terms, variance should be reduced
  var1 <- apply(fixed1, 1, var)
  var2 <- apply(fixed2, 1, var)
  expect_true(all(var2 < var1))
})


# Test tmp.pred() ----

test_that("tmp.pred returns correct dimensions", {
  skip_on_cran()
  # Create mock data and coefficients
  data <- data.frame(
    x1 = c(1, 0, 1),
    x2 = c(0, 1, 0),
    group = c("A", "B", "A")
  )
  coefs <- array(
    rnorm(200),
    dim = c(100, 2, 2),
    dimnames = list(1:100, c("x1", "x2"), c("A", "B"))
  )
  result <- merTools:::tmp.pred(data, coefs, "group")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 100)
})

test_that("tmp.pred warns for new levels", {
  skip_on_cran()
  data <- data.frame(
    x1 = c(1, 0, 1),
    group = c("A", "C", "A")  # "C" is a new level
  )
  coefs <- array(
    rnorm(200),
    dim = c(100, 1, 2),
    dimnames = list(1:100, "x1", c("A", "B"))
  )
  expect_warning(
    merTools:::tmp.pred(data, coefs, "group"),
    "not in the model data"
  )
})

test_that("tmp.pred returns zeros for new levels", {
  skip_on_cran()
  data <- data.frame(
    x1 = c(1, 1),
    group = c("A", "C")  # "C" is a new level
  )
  coefs <- array(
    rep(5, 200),  # Non-zero coefficients
    dim = c(100, 1, 2),
    dimnames = list(1:100, "x1", c("A", "B"))
  )
  result <- suppressWarnings(merTools:::tmp.pred(data, coefs, "group"))
  # Row 1 (level A) should have non-zero predictions
  expect_true(all(result[1, ] != 0))
  # Row 2 (new level C) should have zero predictions
  expect_true(all(result[2, ] == 0))
})


# Test simulate_random_effects() ----

test_that("simulate_random_effects returns named list", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  re_list <- merTools:::simulate_random_effects(
    m1,
    newdata = sleepstudy[1:10, ],
    n.sims = 100
  )
  expect_type(re_list, "list")
  expect_named(re_list, "Subject")
  expect_equal(nrow(re_list$Subject), 10)
  expect_equal(ncol(re_list$Subject), 100)
})

test_that("simulate_random_effects handles multiple grouping factors", {
  skip_on_cran()
  skip_on_ci()
  moddf <- InstEval[sample(rownames(InstEval), 5000), ]
  g1 <- lmer(y ~ lectage + studage + (1 | d) + (1 | s), data = moddf)
  re_list <- merTools:::simulate_random_effects(
    g1,
    newdata = moddf[1:10, ],
    n.sims = 50
  )
  expect_type(re_list, "list")
  expect_true(all(c("d", "s") %in% names(re_list)))
})


# Test combine_components() ----

test_that("combine_components returns correct structure for each which option", {
  skip_on_cran()
  # Create mock data
  fixed_mat <- matrix(rnorm(300), nrow = 3, ncol = 100)
  random_list <- list(
    Subject = matrix(rnorm(300), nrow = 3, ncol = 100)
  )
  sigma_vec <- rep(1, 100)

  # Test which = "full"
  result_full <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = FALSE, which = "full"
  )
  expect_true(is.matrix(result_full))
  expect_equal(dim(result_full), c(3, 100))

  # Test which = "fixed"
  result_fixed <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = FALSE, which = "fixed"
  )
  expect_true(is.matrix(result_fixed))
  expect_equal(dim(result_fixed), c(3, 100))

  # Test which = "random"
  result_random <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = FALSE, which = "random"
  )
  expect_true(is.matrix(result_random))
  expect_equal(dim(result_random), c(3, 100))

  # Test which = "all"
  result_all <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = FALSE, which = "all"
  )
  expect_type(result_all, "list")
  expect_named(result_all, c("yhat", "components"))
  expect_true(is.matrix(result_all$yhat))
})

test_that("combine_components applies residual variance", {
  skip_on_cran()
  set.seed(123)
  fixed_mat <- matrix(rep(0, 300), nrow = 3, ncol = 100)
  random_list <- list(
    Subject = matrix(rep(0, 300), nrow = 3, ncol = 100)
  )
  sigma_vec <- rep(10, 100)  # Large sigma

  result_with <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = TRUE, which = "full"
  )
  result_without <- merTools:::combine_components(
    fixed_mat, random_list, sigma_vec,
    include.resid.var = FALSE, which = "full"
  )

  # With residual variance, results should have larger spread
  expect_gt(sd(result_with), sd(result_without))
})


# Test summarise_predictions() ----

test_that("summarise_predictions returns correct columns", {
  skip_on_cran()
  yhat_arr <- matrix(rnorm(300), nrow = 3, ncol = 100)
  result <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "median",
    predict.type = "linear.prediction",
    N = 3
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_named(result, c("fit", "upr", "lwr"))
})

test_that("summarise_predictions respects stat.type", {
  skip_on_cran()
  set.seed(123)
  yhat_arr <- matrix(rnorm(3000, mean = 100), nrow = 3, ncol = 1000)
  result_median <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "median",
    predict.type = "linear.prediction",
    N = 3
  )
  result_mean <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "mean",
    predict.type = "linear.prediction",
    N = 3
  )
  # Mean and median should be close but not identical
  expect_equal(result_median$fit, result_mean$fit, tolerance = 1)
  expect_false(identical(result_median$fit, result_mean$fit))
})

test_that("summarise_predictions respects level parameter", {
  skip_on_cran()
  set.seed(123)
  yhat_arr <- matrix(rnorm(3000), nrow = 3, ncol = 1000)
  result_80 <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "median",
    predict.type = "linear.prediction",
    N = 3
  )
  result_95 <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.95,
    stat.type = "median",
    predict.type = "linear.prediction",
    N = 3
  )
  # 95% interval should be wider than 80%
  width_80 <- result_80$upr - result_80$lwr
  width_95 <- result_95$upr - result_95$lwr
  expect_true(all(width_95 > width_80))
})

test_that("summarise_predictions applies link function for probability", {
  skip_on_cran()
  d1 <- cbpp
  d1$y <- d1$incidence / d1$size
  gm1 <- glmer(y ~ period + (1 | herd), family = binomial, data = d1, weights = d1$size)

  yhat_arr <- matrix(rnorm(300, mean = 0, sd = 1), nrow = 3, ncol = 100)
  result <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "median",
    predict.type = "probability",
    N = 3,
    merMod = gm1
  )
  # All values should be between 0 and 1 for probability
  expect_true(all(result$fit >= 0 & result$fit <= 1))
  expect_true(all(result$upr >= 0 & result$upr <= 1))
  expect_true(all(result$lwr >= 0 & result$lwr <= 1))
})

test_that("summarise_predictions handles which='all' correctly", {
  skip_on_cran()
  yhat_arr <- matrix(rnorm(300), nrow = 3, ncol = 100)
  pi.comps <- list(
    Subject = matrix(rnorm(300), nrow = 3, ncol = 100),
    fixed = matrix(rnorm(300), nrow = 3, ncol = 100)
  )
  result <- merTools:::summarise_predictions(
    yhat_arr,
    level = 0.8,
    stat.type = "median",
    predict.type = "linear.prediction",
    N = 3,
    which.eff = "all",
    pi.comps = pi.comps
  )
  # Should have effect and obs columns
  expect_true("effect" %in% names(result))
  expect_true("obs" %in% names(result))
  # Should have combined + 2 components = 9 rows (3 obs each)
  expect_equal(nrow(result), 9)
})
