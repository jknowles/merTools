# -----------------------------------------------------------------------------
# Layer 2 numeric regression tests for predictInterval().
#
# These tests pin the *exact* predictInterval() output for a canonical set of
# LMM and GLMM inputs via testthat snapshots. Any change to simulation
# internals that alters these outputs produces a reviewable snapshot diff,
# not a mysterious random-tolerance failure.
#
# Scope: Linux only. Apple-silicon and Windows produce bit-level drift from
# different BLAS/LAPACK, which is benign but would pollute the snapshots.
# Cross-platform numerical equivalence is checked manually by the harness at
# tests/comparisons/predictInterval-regression.R. Linux is the canonical
# platform for the committed snapshot contract.
#
# If a change is intentional, regenerate with:
#   testthat::snapshot_accept("predictInterval-snapshot")
# and inspect the diff carefully before committing.
# -----------------------------------------------------------------------------

SEED  <- 11213
NSIMS <- 1000
TOL   <- 1e-6

snapshot_platform_ok <- function() {
  skip_on_os(c("mac", "windows", "solaris"))
  skip_on_cran()
}

# LMM: random slope + random intercept --------------------------------------

test_that("LMM random-slope + intercept: predictInterval output is stable", {
  snapshot_platform_ok()

  data(sleepstudy, package = "lme4")
  m  <- suppressMessages(
    lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  )
  nd <- sleepstudy[1:10, ]

  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    include.resid.var = FALSE),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    which = "fixed"),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    which = "random"),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    which = "all"),
    style = "json2", tolerance = TOL
  )
})

# LMM: random intercept only ------------------------------------------------

test_that("LMM random-intercept: predictInterval output is stable", {
  snapshot_platform_ok()

  data(sleepstudy, package = "lme4")
  m  <- suppressMessages(
    lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  )
  nd <- sleepstudy[1:10, ]

  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    include.resid.var = FALSE),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    stat = "mean"),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    level = 0.99),
    style = "json2", tolerance = TOL
  )
})

# LMM: edge-case inputs ------------------------------------------------------

test_that("LMM: single-row newdata and variance-adjustment options are stable", {
  snapshot_platform_ok()

  data(sleepstudy, package = "lme4")
  m_slope <- suppressMessages(
    lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  )
  m_int <- suppressMessages(
    lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  )
  nd <- sleepstudy[1:10, ]

  # Single-row newdata
  expect_snapshot_value(
    predictInterval(m_slope, newdata = sleepstudy[1, , drop = FALSE],
                    seed = SEED, n.sims = 500),
    style = "json2", tolerance = TOL
  )

  # ignore.fixed.terms treats a fixed term as known (zero variance)
  expect_snapshot_value(
    predictInterval(m_int, newdata = nd, seed = SEED, n.sims = NSIMS,
                    include.resid.var = FALSE,
                    ignore.fixed.terms = "(Intercept)"),
    style = "json2", tolerance = TOL
  )

  # fix.intercept.variance adjusts intercept variance for RE covariance
  expect_snapshot_value(
    predictInterval(m_int, newdata = nd, seed = SEED, n.sims = NSIMS,
                    fix.intercept.variance = TRUE),
    style = "json2", tolerance = TOL
  )
})

# NOTE: pathological / rank-deficient models (e.g.,
# `Reaction ~ 0 + Days + Days:Subject + (1 | Days)`) are deliberately NOT in
# Layer 2. `mvtnorm::rmvnorm(method = "chol")` on a near-singular covariance
# matrix gives output that differs across BLAS/LAPACK implementations by
# enough to blow past any reasonable snapshot tolerance, even on the same
# OS across R versions. Layer 2 contracts require numerical stability across
# every target runner; if the model cannot deliver that, snapshot it is
# misplaced there. Such models are covered by:
#   - Layer 1 structural invariants in test-predict.R (dim, type, ordering)
#   - Layer 3 manual regression harness at
#     tests/comparisons/predictInterval-regression.R
#
# Rule of thumb: before adding a model to this file, ask whether the design
# matrix is rank-deficient or near-singular. If yes, skip Layer 2.

# GLMM: binomial ------------------------------------------------------------

test_that("GLMM binomial: predictInterval output is stable", {
  snapshot_platform_ok()

  data(cbpp, package = "lme4")
  m <- suppressMessages(
    glmer(
      cbind(incidence, size - incidence) ~ period + (1 | herd),
      data = cbpp, family = binomial
    )
  )
  nd <- cbpp[1:5, ]

  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    include.resid.var = FALSE, type = "linear.prediction"),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    predictInterval(m, newdata = nd, seed = SEED, n.sims = NSIMS,
                    include.resid.var = FALSE, type = "probability"),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    suppressWarnings(predictInterval(
      m, newdata = nd, seed = SEED, n.sims = NSIMS,
      include.resid.var = TRUE, type = "linear.prediction"
    )),
    style = "json2", tolerance = TOL
  )
  expect_snapshot_value(
    suppressWarnings(predictInterval(
      m, newdata = nd, seed = SEED, n.sims = NSIMS,
      include.resid.var = TRUE, type = "probability"
    )),
    style = "json2", tolerance = TOL
  )
})
