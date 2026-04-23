# -----------------------------------------------------------------------------
# Layer 1 structural invariants for predictInterval().
#
# These are properties that MUST hold for every valid predictInterval() output,
# regardless of seed, simulation noise, platform, model type, or option
# combination. If any of these fail, the function is broken — not the
# tolerance and not the test.
#
# These tests run on every CI platform (no skip_on_os). They complement:
#
#   - test-predictInterval-snapshot.R  (Layer 2, Linux-only, pinned numerics)
#   - tests/comparisons/*              (Layer 3, manual regression harness)
#
# Adding a new invariant: derive the property from the function's contract
# (documentation, mathematical definition, or design intent). Do NOT snapshot
# numeric output here; that belongs in Layer 2. Invariants should be true by
# definition — e.g., "upr >= fit" is true by quantile ordering.
# -----------------------------------------------------------------------------

SEED <- 11213

# Build a small, well-conditioned fixture set covering the main algorithmic
# branches: random-slope LMM, random-intercept LMM, binomial GLMM.
make_fixtures <- function() {
  data(sleepstudy, package = "lme4")
  data(cbpp,       package = "lme4")

  list(
    lmm_slope = list(
      fit = suppressMessages(
        lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
      ),
      newdata = sleepstudy[1:20, ]
    ),
    lmm_int = list(
      fit = suppressMessages(
        lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
      ),
      newdata = sleepstudy[1:20, ]
    ),
    glmm_bin = list(
      fit = suppressMessages(
        glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
      ),
      newdata = cbpp[1:10, ]
    )
  )
}

# --------------------------------------------------------------------------
# Invariant 1: output is a data.frame with the expected column names.
# --------------------------------------------------------------------------

test_that("output is a data.frame with columns fit, upr, lwr", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    out <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 200)
    )
    expect_true(is.data.frame(out), info = nm)
    expect_identical(names(out), c("fit", "upr", "lwr"), info = nm)
  }
})

# --------------------------------------------------------------------------
# Invariant 2: interval ordering — upr >= fit >= lwr for every row.
# Ties are acceptable; strict inequality is not guaranteed (e.g., single-
# level factor with zero variance contribution).
# --------------------------------------------------------------------------

test_that("upr >= fit >= lwr for every row", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    out <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 200)
    )
    expect_true(all(out$upr >= out$fit), info = nm)
    expect_true(all(out$fit >= out$lwr), info = nm)
  }
})

# --------------------------------------------------------------------------
# Invariant 3: row-wise correspondence with newdata for which != "all".
# --------------------------------------------------------------------------

test_that("nrow(output) matches nrow(newdata) for which = full/fixed/random", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    nd <- fx[[nm]]$newdata
    for (w in c("full", "fixed", "random")) {
      out <- suppressWarnings(
        predictInterval(fx[[nm]]$fit, newdata = nd,
                        seed = SEED, n.sims = 200, which = w)
      )
      expect_equal(nrow(out), nrow(nd),
                   info = paste0(nm, " / which = ", w))
    }
  }
})

# --------------------------------------------------------------------------
# Invariant 4: no NAs, NaNs, or Infs in output for well-conditioned inputs.
# --------------------------------------------------------------------------

test_that("no NA / NaN / Inf in output for well-conditioned inputs", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    out <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 200)
    )
    expect_true(all(is.finite(out$fit)), info = nm)
    expect_true(all(is.finite(out$upr)), info = nm)
    expect_true(all(is.finite(out$lwr)), info = nm)
  }
})

# --------------------------------------------------------------------------
# Invariant 5: seed reproducibility — same seed produces bit-identical
# output, different seeds produce different output.
# --------------------------------------------------------------------------

test_that("same seed produces identical output; different seeds do not", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    a <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 200)
    )
    b <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 200)
    )
    c <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED + 1L, n.sims = 200)
    )
    expect_identical(a, b, info = paste0(nm, " / same seed"))
    expect_false(identical(a, c), info = paste0(nm, " / different seed"))
  }
})

# --------------------------------------------------------------------------
# Invariant 6: quantile monotonicity — a higher `level` produces intervals
# that contain the intervals at any lower `level`, since both are computed
# from the same simulation draws.
# --------------------------------------------------------------------------

test_that("higher level produces wider-or-equal intervals on same draws", {
  skip_on_cran()
  fx <- make_fixtures()
  for (nm in names(fx)) {
    i80 <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 500, level = 0.80)
    )
    i99 <- suppressWarnings(
      predictInterval(fx[[nm]]$fit, newdata = fx[[nm]]$newdata,
                      seed = SEED, n.sims = 500, level = 0.99)
    )
    expect_true(all(i99$upr >= i80$upr), info = paste0(nm, " / upr"))
    expect_true(all(i99$lwr <= i80$lwr), info = paste0(nm, " / lwr"))
  }
})

# --------------------------------------------------------------------------
# Invariant 7: GLMM probability predictions must lie in [0, 1].
# --------------------------------------------------------------------------

test_that("GLMM type='probability' output is in [0, 1]", {
  skip_on_cran()
  fx <- make_fixtures()
  out <- suppressWarnings(
    predictInterval(fx$glmm_bin$fit, newdata = fx$glmm_bin$newdata,
                    seed = SEED, n.sims = 500,
                    type = "probability", include.resid.var = FALSE)
  )
  expect_true(all(out$fit >= 0 & out$fit <= 1))
  expect_true(all(out$upr >= 0 & out$upr <= 1))
  expect_true(all(out$lwr >= 0 & out$lwr <= 1))
})

# --------------------------------------------------------------------------
# Regression test: GitHub issue #101
# predictInterval() must not crash on no-intercept random-slope models.
# --------------------------------------------------------------------------

test_that("no-intercept random-slope models do not crash (#101)", {
  skip_on_cran()
  data(sleepstudy, package = "lme4")
  set.seed(SEED)
  sleepstudy$Test <- rep(
    sample(c(TRUE, FALSE), length(unique(sleepstudy$Subject)), replace = TRUE),
    each = 10
  )
  fit <- suppressMessages(
    lmer(Reaction ~ Days:Test + (0 + Days | Subject), data = sleepstudy)
  )
  out <- suppressMessages(
    predictInterval(fit, newdata = head(sleepstudy), n.sims = 200, seed = SEED)
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 6L)
  expect_named(out, c("fit", "upr", "lwr"))
  expect_true(all(out$upr >= out$fit))
  expect_true(all(out$fit >= out$lwr))
})
