# ---------------------------------------------------------------------------
# predictInterval-regression.R
#
# Cross-version numeric regression harness for predictInterval().
#
# Generates a deterministic bundle of predictInterval() outputs for a pinned
# set of canonical models, inputs, and argument combinations. Two package
# versions can then be diffed bit-for-bit, making it easy to confirm that a
# refactor did not change numeric behavior for users who rely on a fixed seed.
#
# This script is NOT part of R CMD check. It is a manual tool invoked
# explicitly when someone touches simulation internals.
#
# Usage:
#
#   # (1) Generate an output bundle for one package version
#   Rscript tests/comparisons/predictInterval-regression.R harness \
#           <pkg_path> <output_rds>
#
#   # (2) Diff two bundles
#   Rscript tests/comparisons/predictInterval-regression.R diff \
#           <rds_a> <rds_b>
#
# Typical workflow (current checkout vs. origin/master):
#
#   git worktree add /tmp/mT-old origin/master
#   Rscript tests/comparisons/predictInterval-regression.R harness \
#           /tmp/mT-old  /tmp/old.rds
#   Rscript tests/comparisons/predictInterval-regression.R harness \
#           .           /tmp/new.rds
#   Rscript tests/comparisons/predictInterval-regression.R diff  \
#           /tmp/old.rds /tmp/new.rds
#   git worktree remove /tmp/mT-old
#
# Interpretation:
#   - All cases IDENTICAL (max_abs = 0)  -> refactor preserves behavior.
#   - Numeric differences                -> behavior changed; inspect whether
#                                           the change is intentional
#                                           (e.g., the binomial GLMM residual
#                                           simulation fix introduced in
#                                           merTools 0.9.0 legitimately alters
#                                           the two glmm_bin_with_resid_*
#                                           cases).
#
# Requires: pkgload, lme4.
# ---------------------------------------------------------------------------

SEED  <- 11213
NSIMS <- 1000

# -- Harness -----------------------------------------------------------------

run_harness <- function(pkg_path, out_path) {
  # Pin RNG so R version does not perturb stream behavior
  RNGversion("4.1.0")
  RNGkind(kind = "Mersenne-Twister",
          normal.kind = "Inversion",
          sample.kind = "Rejection")

  suppressPackageStartupMessages({
    library(pkgload)
    load_all(pkg_path, quiet = TRUE)
    library(lme4)
  })

  data(sleepstudy, package = "lme4")
  data(cbpp,       package = "lme4")

  cases <- list()

  # LMM: random slope + random intercept
  m_slope <- suppressMessages(
    lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  )
  cases$lmm_slope_default <- predictInterval(
    m_slope, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS)
  cases$lmm_slope_no_resid <- predictInterval(
    m_slope, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, include.resid.var = FALSE)
  cases$lmm_slope_fixed_only <- predictInterval(
    m_slope, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, which = "fixed")
  cases$lmm_slope_random_only <- predictInterval(
    m_slope, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, which = "random")
  cases$lmm_slope_all <- predictInterval(
    m_slope, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, which = "all")

  # LMM: random intercept only
  m_int <- suppressMessages(
    lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
  )
  cases$lmm_int_default <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS)
  cases$lmm_int_no_resid <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, include.resid.var = FALSE)
  cases$lmm_int_mean_stat <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, stat = "mean")
  cases$lmm_int_level_99 <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, level = 0.99)

  # LMM: single-row newdata
  cases$lmm_single_row <- predictInterval(
    m_slope, newdata = sleepstudy[1, , drop = FALSE],
    seed = SEED, n.sims = 500)

  # LMM: ignore.fixed.terms + fix.intercept.variance
  cases$lmm_ignore_fixed <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, include.resid.var = FALSE,
    ignore.fixed.terms = "(Intercept)")
  cases$lmm_fix_intercept_var <- predictInterval(
    m_int, newdata = sleepstudy[1:10, ],
    seed = SEED, n.sims = NSIMS, fix.intercept.variance = TRUE)

  # GLMM: binomial (the two with_resid cases differ intentionally in 0.9.0+
  # because of the n-trial binomial residual simulation fix)
  gm_bin <- suppressMessages(
    glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
          data = cbpp, family = binomial)
  )
  cases$glmm_bin_no_resid_linear <- predictInterval(
    gm_bin, newdata = cbpp[1:5, ],
    seed = SEED, n.sims = NSIMS,
    include.resid.var = FALSE, type = "linear.prediction")
  cases$glmm_bin_no_resid_prob <- predictInterval(
    gm_bin, newdata = cbpp[1:5, ],
    seed = SEED, n.sims = NSIMS,
    include.resid.var = FALSE, type = "probability")
  cases$glmm_bin_with_resid_linear <- suppressWarnings(predictInterval(
    gm_bin, newdata = cbpp[1:5, ],
    seed = SEED, n.sims = NSIMS,
    include.resid.var = TRUE, type = "linear.prediction"))
  cases$glmm_bin_with_resid_prob <- suppressWarnings(predictInterval(
    gm_bin, newdata = cbpp[1:5, ],
    seed = SEED, n.sims = NSIMS,
    include.resid.var = TRUE, type = "probability"))

  attr(cases, "metadata") <- list(
    pkg_path        = normalizePath(pkg_path),
    package_version = as.character(utils::packageVersion("merTools")),
    R_version       = R.version.string,
    sysname         = Sys.info()[["sysname"]],
    machine         = Sys.info()[["machine"]],
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed            = SEED,
    n_sims          = NSIMS,
    n_cases         = length(cases)
  )

  saveRDS(cases, out_path)
  cat("Saved", length(cases), "cases to", out_path, "\n")
  cat("Package version:", attr(cases, "metadata")$package_version, "\n")
  invisible(cases)
}

# -- Diff --------------------------------------------------------------------

summarize_case <- function(ma, mb) {
  cols <- intersect(names(ma), names(mb))
  if (length(cols) == 0L) {
    return(data.frame(col = "<no common cols>",
                      max_abs = NA_real_, rel_err = NA_real_))
  }
  do.call(rbind, lapply(cols, function(col) {
    x <- suppressWarnings(as.numeric(ma[[col]]))
    y <- suppressWarnings(as.numeric(mb[[col]]))
    if (length(x) != length(y)) {
      return(data.frame(col = col, max_abs = NA_real_, rel_err = NA_real_))
    }
    d <- abs(x - y)
    s <- pmax(abs(x), abs(y), 1e-12)
    data.frame(col = col,
               max_abs = max(d, na.rm = TRUE),
               rel_err = max(d / s, na.rm = TRUE))
  }))
}

classify <- function(max_abs) {
  if (!is.finite(max_abs))  return("<all NA>")
  if (max_abs < 1e-10)      return("IDENTICAL")
  if (max_abs < 1e-6)       return("numeric-equal")
  if (max_abs < 1e-3)       return("tiny drift")
  if (max_abs < 0.5)        return("small drift")
  if (max_abs < 5)          return("moderate drift")
                            return("LARGE DIFF")
}

run_diff <- function(rds_a, rds_b) {
  a <- readRDS(rds_a); b <- readRDS(rds_b)
  meta_a <- attr(a, "metadata"); meta_b <- attr(b, "metadata")

  cat(sprintf("A: %s  (%s)\n", rds_a, meta_a$package_version %||% "?"))
  cat(sprintf("B: %s  (%s)\n", rds_b, meta_b$package_version %||% "?"))
  cat(sprintf("R: %s\n\n", meta_a$R_version %||% R.version.string))

  worst <- 0
  for (nm in union(names(a), names(b))) {
    if (is.null(a[[nm]]) || is.null(b[[nm]])) {
      cat(sprintf("%-32s [missing in one bundle]\n", nm)); next
    }
    ma <- a[[nm]]; mb <- b[[nm]]
    if (!identical(dim(ma), dim(mb))) {
      cat(sprintf("%-32s [structure differs]\n", nm)); next
    }
    r <- summarize_case(ma, mb)
    maxa <- max(r$max_abs, na.rm = TRUE)
    maxr <- max(r$rel_err, na.rm = TRUE)
    worst <- max(worst, if (is.finite(maxa)) maxa else 0)
    cat(sprintf("%-32s  max_abs=%.3e  max_rel=%.3e  %s\n",
                nm, maxa, maxr, classify(maxa)))
  }
  cat(sprintf("\nOverall max absolute diff: %.6e\n", worst))
  invisible(worst)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# -- Dispatch ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript predictInterval-regression.R harness <pkg> <out.rds>\n",
       "   or: Rscript predictInterval-regression.R diff    <a.rds> <b.rds>")
}
mode <- args[[1]]
if (mode == "harness") {
  if (length(args) != 3L) stop("harness mode requires <pkg_path> <output_rds>")
  run_harness(args[[2]], args[[3]])
} else if (mode == "diff") {
  if (length(args) != 3L) stop("diff mode requires <rds_a> <rds_b>")
  run_diff(args[[2]], args[[3]])
} else {
  stop("Unknown mode: ", mode)
}
