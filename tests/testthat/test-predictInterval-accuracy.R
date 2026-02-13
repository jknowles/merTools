test_that("component reproducibility with seed", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 42
  n.sims <- 100
  
  result1 <- predictInterval(
    m,
    newdata = sleepstudy[1:5, ],
    n.sims = n.sims,
    seed = seed
  )
  
  result2 <- predictInterval(
    m,
    newdata = sleepstudy[1:5, ],
    n.sims = n.sims,
    seed = seed
  )
  
  expect_equal(result1, result2)
})

test_that("single observation predictions (n=1)", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 123
  result <- predictInterval(
    m,
    newdata = sleepstudy[1, ],
    n.sims = 500,
    seed = seed
  )
  
  expect_equal(nrow(result), 1)
  expect_identical(names(result), c("fit", "upr", "lwr"))
  expect_true(result$upr > result$fit)
  expect_true(result$fit > result$lwr)
})

test_that("extreme confidence levels (0.8 vs 0.99)", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 456
  
  result_80 <- predictInterval(
    m,
    newdata = sleepstudy[1:3, ],
    level = 0.8,
    n.sims = 500,
    seed = seed
  )
  
  result_99 <- predictInterval(
    m,
    newdata = sleepstudy[1:3, ],
    level = 0.99,
    n.sims = 500,
    seed = seed
  )
  
  expect_equal(nrow(result_80), 3)
  expect_equal(nrow(result_99), 3)
  
  interval_width_80 <- result_80[,"upr"] - result_80[,"lwr"]
  interval_width_99 <- result_99[,"upr"] - result_99[,"lwr"]
  
  expect_true(all(interval_width_99 > interval_width_80))
})

test_that("GLMM probability predictions bounded [0,1]", {
  data(cbpp, package = "lme4")
  
  m <- glmer(
        cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial
  )
  
  seed <- 789
  result <- predictInterval(
    m,
    newdata = cbpp[1:10, ],
    type = "probability",
    n.sims = 200,
    seed = seed
  )
  
  expect_equal(nrow(result), 10)
  expect_true(all(result[, "fit"] >= 0 & result[, "fit"] <= 1))
  expect_true(all(result[, "lwr"] >= 0 & result[, "lwr"] <= 1))
  expect_true(all(result[, "upr"] >= 0 & result[, "upr"] <= 1))
})

test_that("parallelization does not throw errors", {
  skip_on_cran()
  
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 999
  
  expect_no_error({
    predictInterval(
      m,
      newdata = sleepstudy[1:20, ],
      n.sims = 500,
      seed = seed,
      .parallel = TRUE
    )
  })
})

test_that("all which options return correct dimensions", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 111
  n.sims <- 200
  n.obs <- 5
  
  result_full <- predictInterval(
    m,
    newdata = sleepstudy[1:n.obs, ],
    which = "full",
    n.sims = n.sims,
    seed = seed
  )
  
  result_fixed <- predictInterval(
    m,
    newdata = sleepstudy[1:n.obs, ],
    which = "fixed",
    n.sims = n.sims,
    seed = seed
  )
  
  result_random <- predictInterval(
    m,
    newdata = sleepstudy[1:n.obs, ],
    which = "random",
    n.sims = n.sims,
    seed = seed
  )
  
  result_all <- predictInterval(
    m,
    newdata = sleepstudy[1:n.obs, ],
    which = "all",
    n.sims = n.sims,
    seed = seed
  )
  
  expect_equal(nrow(result_full), n.obs)
  expect_equal(nrow(result_fixed), n.obs)
  expect_equal(nrow(result_random), n.obs)
  
  expect_identical(names(result_full), c("fit", "upr", "lwr"))
  expect_identical(names(result_fixed), c("fit", "upr", "lwr"))
  expect_identical(names(result_random), c("fit", "upr", "lwr"))
  
  expect_true(nrow(result_all) > n.obs)
  expect_true("effect" %in% names(result_all))
})

test_that("large simulation count (n.sims = 10000)", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 222
  n.sims <- 10000
  
  result <- predictInterval(
    m,
    newdata = sleepstudy[1:10, ],
    n.sims = n.sims,
    seed = seed,
    returnSims = TRUE
  )
  
  expect_equal(nrow(result), 10)
  
  attr_sims <- attr(result, "sim.results")
  expect_true(!is.null(attr_sims))
  
  interval_width <- result[, "upr"] - result[, "lwr"]
  expect_true(all(interval_width > 0))
})

test_that("component-by-component consistency with predictInterval wrapper", {
  data(sleepstudy, package = "lme4")
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
  
  seed <- 12345
  n.sims <- 200
  newdata <- sleepstudy[1:5, ]
  
  set.seed(seed)
  full_result <- predictInterval(
    m,
    newdata = newdata,
    n.sims = n.sims,
    seed = seed
  )
  
  sigma_vec <- merTools:::simulate_residual_variance(m, n.sims)
  random_list <- merTools:::simulate_random_effects(m, newdata, n.sims, seed = seed)
  fixed_mat <- merTools:::simulate_fixed_effects(m, newdata, n.sims, seed = seed)
  
  combined <- merTools:::combine_components(
    fixed_mat = fixed_mat,
    random_list = random_list,
    sigma_vec = sigma_vec,
    include.resid.var = TRUE,
    which = "full"
  )
  
  combined_df <- as.data.frame(combined)
  names(combined_df) <- c("fit", "upr", "lwr")
  
  expect_equal(full_result$fit, combined_df$fit, tolerance = sd(full_result$fit)/100)
  expect_equal(full_result$upr, combined_df$upr, tolerance = sd(full_result$upr)/100)
  expect_equal(full_result$lwr, combined_df$lwr, tolerance = sd(full_result$lwr)/100)
})
