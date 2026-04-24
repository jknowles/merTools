# -----------------------------------------------------------------------------
# Regression tests for issue #118:
# predictInterval() fails on models with multiple random-effect term blocks
# for the same grouping factor (double-bar syntax, explicit splits, mixed
# correlated + uncorrelated). Root cause: lme4 returns `postVar` as a list
# of arrays in those cases, and the indexing code assumed a 3-D array.
# Fixed by combine_postvar_blocks() which normalizes to a block-diagonal
# array before the level-filtering and rmvnorm path runs.
#
# Structural invariants only — numeric output for uncorrelated specs is
# sensitive to lme4's internal parameterization and is not pinned here.
# -----------------------------------------------------------------------------

SEED <- 11213
N_SIMS <- 500

make_sleep_fixture <- function() {
  data(sleepstudy, package = "lme4")
  set.seed(SEED)
  sleepstudy$av <- rnorm(nrow(sleepstudy))
  sleepstudy
}

expect_valid_predictInterval <- function(pi, n_rows) {
  expect_s3_class(pi, "data.frame")
  expect_equal(nrow(pi), n_rows)
  expect_named(pi, c("fit", "upr", "lwr"))
  expect_true(all(is.finite(unlist(pi))))
  expect_true(all(pi$upr >= pi$fit))
  expect_true(all(pi$fit >= pi$lwr))
}

# 1. Exact reproducer from issue #118 ------------------------------------------
test_that("predictInterval works with double-bar uncorrelated random slopes (#118)", {
  skip_on_cran()
  sleep <- make_sleep_fixture()
  fm <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (Days + av || Subject), data = sleep)
  ))
  pi <- suppressWarnings(
    predictInterval(fm, newdata = sleep, n.sims = N_SIMS, seed = SEED)
  )
  expect_valid_predictInterval(pi, nrow(sleep))
})

# 2. Double-bar with intercept -------------------------------------------------
test_that("predictInterval works with (1 + Days || Subject)", {
  skip_on_cran()
  data(sleepstudy, package = "lme4")
  fm <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (1 + Days || Subject), data = sleepstudy)
  ))
  pi <- suppressWarnings(
    predictInterval(fm, newdata = sleepstudy, n.sims = N_SIMS, seed = SEED)
  )
  expect_valid_predictInterval(pi, nrow(sleepstudy))
})

# 3. Explicit split equivalent -------------------------------------------------
test_that("predictInterval works with explicit split (1|Subject) + (0 + Days|Subject)", {
  skip_on_cran()
  data(sleepstudy, package = "lme4")
  fm <- suppressMessages(suppressWarnings(
    lme4::lmer(
      Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
      data = sleepstudy
    )
  ))
  pi <- suppressWarnings(
    predictInterval(fm, newdata = sleepstudy, n.sims = N_SIMS, seed = SEED)
  )
  expect_valid_predictInterval(pi, nrow(sleepstudy))
})

# 4. Mixed correlated + uncorrelated -------------------------------------------
test_that("predictInterval works with mixed (Days|Subject) + (0 + av|Subject)", {
  skip_on_cran()
  sleep <- make_sleep_fixture()
  fm <- suppressMessages(suppressWarnings(
    lme4::lmer(
      Reaction ~ Days + (Days | Subject) + (0 + av | Subject),
      data = sleep
    )
  ))
  pi <- suppressWarnings(
    predictInterval(fm, newdata = sleep, n.sims = N_SIMS, seed = SEED)
  )
  expect_valid_predictInterval(pi, nrow(sleep))
})

# 5. GLMM with uncorrelated random slopes --------------------------------------
test_that("predictInterval works on GLMM with uncorrelated random slopes", {
  skip_on_cran()
  data(cbpp, package = "lme4")
  cbpp$period_num <- as.numeric(cbpp$period)
  gm <- suppressMessages(suppressWarnings(
    lme4::glmer(
      cbind(incidence, size - incidence) ~ period + (1 + period_num || herd),
      data = cbpp, family = binomial
    )
  ))
  pi <- suppressWarnings(
    predictInterval(
      gm, newdata = cbpp, n.sims = N_SIMS, seed = SEED,
      type = "linear.prediction", include.resid.var = FALSE
    )
  )
  expect_valid_predictInterval(pi, nrow(cbpp))
})

# 6. Correlated-only case (regression guard) -----------------------------------
test_that("correlated-only random slopes still work (regression guard)", {
  skip_on_cran()
  data(sleepstudy, package = "lme4")
  fm <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  ))
  pi <- suppressWarnings(
    predictInterval(fm, newdata = sleepstudy, n.sims = N_SIMS, seed = SEED)
  )
  expect_valid_predictInterval(pi, nrow(sleepstudy))
})

# 7. Unit test on combine_postvar_blocks() -------------------------------------
test_that("combine_postvar_blocks builds correct block-diagonal array", {
  # Two blocks: a 2x2 with level-specific covariance, and a 1x1 scalar.
  # Three levels. Fill with distinctive values so block-diagonal structure
  # is verifiable element-wise.
  n_levels <- 3L
  block1 <- array(0, dim = c(2L, 2L, n_levels))
  block2 <- array(0, dim = c(1L, 1L, n_levels))
  for (k in seq_len(n_levels)) {
    block1[, , k] <- matrix(c(k, 0.1 * k, 0.1 * k, 2 * k), 2L, 2L)
    block2[, , k] <- matrix(10 * k, 1L, 1L)
  }

  # Mock reMeans: 3 levels x 3 columns, names carried to dimnames of output.
  reMeans <- matrix(0, nrow = n_levels, ncol = 3L,
                    dimnames = list(c("L1", "L2", "L3"),
                                    c("(Intercept)", "x", "y")))

  out <- combine_postvar_blocks(list(block1, block2), reMeans)

  expect_true(is.array(out))
  expect_equal(dim(out), c(3L, 3L, 3L))
  expect_equal(dimnames(out)[[1]], c("(Intercept)", "x", "y"))
  expect_equal(dimnames(out)[[3]], c("L1", "L2", "L3"))

  for (k in seq_len(n_levels)) {
    # Block 1 occupies rows/cols 1:2
    expect_equal(out[1:2, 1:2, k], block1[, , k], ignore_attr = TRUE)
    # Block 2 occupies row/col 3
    expect_equal(out[3, 3, k], 10 * k, ignore_attr = TRUE)
    # Off-diagonal entries between blocks are zero
    expect_equal(out[1:2, 3, k], c(0, 0), ignore_attr = TRUE)
    expect_equal(out[3, 1:2, k], c(0, 0), ignore_attr = TRUE)
  }
})

# 8. Helper raises informative error on malformed input ------------------------
test_that("combine_postvar_blocks stops when block dims don't sum to ncol(reMeans)", {
  block <- array(0, dim = c(1L, 1L, 2L))
  reMeans <- matrix(0, nrow = 2L, ncol = 3L,
                    dimnames = list(c("L1", "L2"), c("a", "b", "c")))
  expect_error(
    combine_postvar_blocks(list(block), reMeans),
    "Unexpected lme4 postVar structure"
  )
})
