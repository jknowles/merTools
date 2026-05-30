# ----------------------------------------------------------------------------
# brms vs. lme4 + merTools::predictInterval(): prediction intervals compared
# ----------------------------------------------------------------------------
#
# Since merTools was written, brms has become a de-facto standard for fully
# Bayesian multilevel modeling. This script puts the fast, simulation-based
# intervals from `predictInterval()` next to the "gold standard" posterior
# intervals from brms across three studies:
#
#   STUDY 1  Gaussian random-slopes model (hsb: math achievement)
#            -- point estimates, interval agreement, out-of-sample coverage,
#               speed, and behaviour for entirely new groups (unseen schools).
#   STUDY 2  The same model WITHOUT the contextual fixed effect (meanses)
#            -- shows that predictInterval's new-group intervals get too narrow
#               (and the gap to brms widens) once a fixed effect no longer
#               absorbs the between-school variance.
#   STUDY 3  Binomial GLMM (grouseticks) compared on the probability scale
#            -- predictInterval(type = "probability") vs brms posterior_epred,
#               plus a calibration check.
#
# Requires: merTools (>= 1.0.0), lme4, brms, ggplot2.
#
# Output and the (cached) brms fits go to `MERTOOLS_BENCH_DIR` if set, else a
# per-session temp directory. Each brms fit is cached, so only the first run
# pays the Stan compile + sampling cost.
# ----------------------------------------------------------------------------

if (!requireNamespace("brms", quietly = TRUE)) {
  stop("This example needs the 'brms' package. install.packages('brms')")
}
suppressPackageStartupMessages({
  library(merTools)
  library(lme4)
  library(brms)
  library(ggplot2)
})

SEED       <- 11213
N_SIMS     <- 1000
NOMINAL    <- c(0.50, 0.80, 0.90, 0.95)
N_NEW_GRPS <- 6
TEST_FRAC  <- 0.20
BRM_ITER   <- 4000
BRM_WARMUP <- 1000

outdir <- Sys.getenv("MERTOOLS_BENCH_DIR", unset = "")
if (!nzchar(outdir)) outdir <- file.path(tempdir(), "brms_vs_predictInterval")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", outdir)
set.seed(SEED)

## ===========================================================================
## Shared helpers
## ===========================================================================

# Standardised brms fit with on-disk caching.
fit_brm <- function(form, data, family, cache) {
  brm(form, data = data, family = family, chains = 4, iter = BRM_ITER,
      warmup = BRM_WARMUP, seed = SEED, refresh = 0, silent = 2,
      control = list(adapt_delta = 0.95),
      file = file.path(outdir, cache), file_refit = "on_change")
}

# (N_SIMS x nrow(newdata)) matrix of predictive draws from predictInterval(),
# oriented like brms::posterior_predict() (draws x obs).
mt_draws <- function(fit, newdata, type = "linear.prediction",
                     include.resid.var = TRUE) {
  pi <- suppressWarnings(predictInterval(
    fit, newdata = newdata, level = 0.95, n.sims = N_SIMS, seed = SEED,
    type = type, include.resid.var = include.resid.var, returnSims = TRUE))
  sims <- attr(pi, "sim.results")
  if (nrow(sims) == nrow(newdata)) sims <- t(sims)
  sims
}

qband <- function(draws, level) {
  a <- (1 - level) / 2
  data.frame(fit = apply(draws, 2, median),
             lwr = apply(draws, 2, function(x) unname(quantile(x, a))),
             upr = apply(draws, 2, function(x) unname(quantile(x, 1 - a))))
}

coverage <- function(draws, y, levels) {
  vapply(levels, function(p) {
    b <- qband(draws, p); mean(y >= b$lwr & y <= b$upr)
  }, numeric(1))
}

# Full report for one Gaussian study. Returns tables + draws (for figures).
report_gaussian <- function(label, m_lme, m_brm, test_seen, test_new) {
  cat("\n##########  ", label, "  ##########\n", sep = "")
  t_mt <- system.time({
    mt_seen <- mt_draws(m_lme, test_seen)
    mt_new  <- mt_draws(m_lme, test_new)
  })
  t_pp <- system.time({
    pp_seen <- posterior_predict(m_brm, newdata = test_seen, allow_new_levels = TRUE)
    pp_new  <- posterior_predict(m_brm, newdata = test_new, allow_new_levels = TRUE)
  })

  # 1. point estimates (conditional means)
  lme_mean <- predict(m_lme, newdata = test_seen)
  brm_mean <- colMeans(posterior_epred(m_brm, newdata = test_seen,
                                       allow_new_levels = TRUE))
  pe <- lme_mean - brm_mean
  cat(sprintf("\n[1] Point estimates  cor=%.4f  mean|diff|=%.4f  RMSE=%.4f  (y SD=%.2f)\n",
              cor(lme_mean, brm_mean), mean(abs(pe)), sqrt(mean(pe^2)),
              sd(test_seen[[all.vars(formula(m_lme))[1]]])))

  # 2. interval agreement (90%)
  b_mt <- qband(mt_seen, 0.90); b_brm <- qband(pp_seen, 0.90)
  cat(sprintf("[2] 90%% PI width  merTools=%.2f brms=%.2f  cor(lwr/upr)=%.3f/%.3f  mean|endpt diff|=%.3f\n",
              median(b_mt$upr - b_mt$lwr), median(b_brm$upr - b_brm$lwr),
              cor(b_mt$lwr, b_brm$lwr), cor(b_mt$upr, b_brm$upr),
              mean(abs(c(b_mt$lwr - b_brm$lwr, b_mt$upr - b_brm$upr)))))

  # 3. coverage (seen) and 4. new groups
  y_seen <- test_seen[[all.vars(formula(m_lme))[1]]]
  y_new  <- test_new[[all.vars(formula(m_lme))[1]]]
  cov_seen <- data.frame(nominal = NOMINAL,
                         merTools = coverage(mt_seen, y_seen, NOMINAL),
                         brms     = coverage(pp_seen, y_seen, NOMINAL))
  cov_new  <- data.frame(nominal = NOMINAL,
                         merTools = coverage(mt_new, y_new, NOMINAL),
                         brms     = coverage(pp_new, y_new, NOMINAL))
  cat("[3] Coverage, held-out students in SEEN schools:\n")
  print(cov_seen, row.names = FALSE, digits = 3)
  sd_school <- attr(lme4::VarCorr(m_lme)[[1]], "stddev")[["(Intercept)"]]
  cat(sprintf("[4] New schools: school-intercept SD=%.2f vs residual SD=%.2f\n",
              sd_school, sigma(m_lme)))
  print(cov_new, row.names = FALSE, digits = 3)

  invisible(list(point = data.frame(merTools = lme_mean, brms = brm_mean),
                 b_mt = b_mt, b_brm = b_brm,
                 cov_seen = cov_seen, cov_new = cov_new,
                 sd_school = sd_school,
                 timing = c(mt = t_mt[["elapsed"]], pp = t_pp[["elapsed"]])))
}

## ===========================================================================
## Data: hsb, with a train / test-seen / test-new split
## ===========================================================================
data(hsb, package = "merTools")
hsb$schid <- factor(hsb$schid)
new_schools <- sample(levels(hsb$schid), N_NEW_GRPS)
hsb$.set <- "train"
hsb$.set[hsb$schid %in% new_schools] <- "test_new"
seen_rows <- which(hsb$.set == "train")
hsb$.set[sample(seen_rows, round(TEST_FRAC * length(seen_rows)))] <- "test_seen"
train     <- droplevels(hsb[hsb$.set == "train", ])
test_seen <- hsb[hsb$.set == "test_seen" & hsb$schid %in% levels(train$schid), ]
test_new  <- hsb[hsb$.set == "test_new", ]
cat(sprintf("\nhsb split -> train %d / test_seen %d / test_new %d (%d new schools)\n",
            nrow(train), nrow(test_seen), nrow(test_new), length(new_schools)))

## ===========================================================================
## STUDY 1 -- Gaussian random slopes, WITH contextual fixed effect (meanses)
## ===========================================================================
f1 <- mathach ~ ses + meanses + (ses | schid)
t_lme1 <- system.time(m_lme1 <- lmer(f1, data = train))
brm1_cached <- file.exists(file.path(outdir, "brms_hsb_meanses.rds"))
t_brm1 <- system.time(m_brm1 <- fit_brm(f1, train, gaussian(), "brms_hsb_meanses"))
res1 <- report_gaussian("STUDY 1: ses + meanses + (ses|schid)",
                        m_lme1, m_brm1, test_seen, test_new)
cat("\n=== STUDY 1 timing (wall-clock) ===\n")
timing <- data.frame(
  step = c("lme4 fit", "predictInterval (seen+new)",
           "brms fit (compile+sample)", "posterior_predict (seen+new)"),
  seconds = c(t_lme1[["elapsed"]], res1$timing[["mt"]],
              t_brm1[["elapsed"]], res1$timing[["pp"]]))
print(timing, row.names = FALSE, digits = 3)
mt_total  <- t_lme1[["elapsed"]] + res1$timing[["mt"]]
brm_total <- t_brm1[["elapsed"]] + res1$timing[["pp"]]
if (brm1_cached) {
  cat("(brms fit was loaded from cache; delete it to time a fresh compile + sample.)\n")
} else {
  cat(sprintf("end-to-end: merTools %.1fs (lme4 + predictInterval) vs brms %.1fs  (~%.0fx faster)\n",
              mt_total, brm_total, brm_total / mt_total))
}

## ===========================================================================
## STUDY 2 -- same model WITHOUT meanses: the new-group contrast
## ===========================================================================
f2 <- mathach ~ ses + (ses | schid)
m_lme2 <- lmer(f2, data = train)
m_brm2 <- fit_brm(f2, train, gaussian(), "brms_hsb_nomeanses")
res2 <- report_gaussian("STUDY 2: ses + (ses|schid)  [no meanses]",
                        m_lme2, m_brm2, test_seen, test_new)

cat("\n=== New-group coverage: contextual effect present vs absent ===\n")
contrast <- data.frame(
  model = c("with meanses", "without meanses"),
  school_SD = c(res1$sd_school, res2$sd_school),
  merTools_90 = c(res1$cov_new$merTools[NOMINAL == 0.90],
                  res2$cov_new$merTools[NOMINAL == 0.90]),
  brms_90 = c(res1$cov_new$brms[NOMINAL == 0.90],
              res2$cov_new$brms[NOMINAL == 0.90]))
contrast$gap_90 <- contrast$brms_90 - contrast$merTools_90
print(contrast, row.names = FALSE, digits = 3)
cat("Dropping meanses raises the school SD; predictInterval (which omits the\n",
    "school effect for unseen schools) then under-covers more, widening the gap\n",
    "to brms (which samples a new school effect).\n", sep = "")

## ===========================================================================
## STUDY 3 -- Binomial GLMM (grouseticks), compared on the probability scale
## ===========================================================================
data(grouseticks, package = "lme4")
grouseticks$TICKS_BIN <- as.integer(grouseticks$TICKS >= 1)
grouseticks$cHEIGHT   <- as.numeric(scale(grouseticks$HEIGHT))
gk <- grouseticks
gk$.set <- ifelse(runif(nrow(gk)) < TEST_FRAC, "test", "train")
gk_tr <- droplevels(gk[gk$.set == "train", ])
gk_te <- gk[gk$.set == "test" & gk$BROOD %in% levels(gk_tr$BROOD)
            & gk$LOCATION %in% unique(gk_tr$LOCATION), ]

gf <- TICKS_BIN ~ YEAR + cHEIGHT + (1 | BROOD) + (1 | LOCATION)
g_lme <- glmer(gf, family = binomial, data = gk_tr,
               control = glmerControl(optimizer = "bobyqa"))
g_brm <- fit_brm(gf, gk_tr, bernoulli(), "brms_grouseticks")

# probability draws (no Bernoulli noise): predictInterval prob vs brms epred
mt_pr <- mt_draws(g_lme, gk_te, type = "probability", include.resid.var = FALSE)
if (max(mt_pr) > 1 || min(mt_pr) < 0) mt_pr <- plogis(mt_pr)  # guard: ensure prob scale
brm_pr <- posterior_epred(g_brm, newdata = gk_te, allow_new_levels = TRUE)

cat("\n##########  STUDY 3: grouseticks binomial GLMM (probability scale)  ##########\n")
mt_p <- apply(mt_pr, 2, median); brm_p <- apply(brm_pr, 2, median)
cat(sprintf("[1] Predicted probabilities  cor=%.4f  mean|diff|=%.4f\n",
            cor(mt_p, brm_p), mean(abs(mt_p - brm_p))))
bp_mt <- qband(mt_pr, 0.90); bp_brm <- qband(brm_pr, 0.90)
cat(sprintf("[2] 90%% interval on P(tick)  cor(lwr/upr)=%.3f/%.3f  mean|endpt diff|=%.4f\n",
            cor(bp_mt$lwr, bp_brm$lwr), cor(bp_mt$upr, bp_brm$upr),
            mean(abs(c(bp_mt$lwr - bp_brm$lwr, bp_mt$upr - bp_brm$upr)))))
# calibration: bin predicted prob, compare to observed frequency
calib <- function(p, y, k = 5) {
  b <- cut(p, quantile(p, seq(0, 1, length.out = k + 1)), include.lowest = TRUE)
  data.frame(bin = seq_len(nlevels(b)),
             pred = tapply(p, b, mean), obs = tapply(y, b, mean),
             n = as.integer(table(b)))
}
cat("[3] Calibration (predicted vs observed P(tick) by quintile):\n")
cal_mt  <- calib(mt_p,  gk_te$TICKS_BIN); cal_mt$method  <- "merTools"
cal_brm <- calib(brm_p, gk_te$TICKS_BIN); cal_brm$method <- "brms"
print(rbind(cal_mt, cal_brm)[, c("method", "bin", "pred", "obs", "n")],
      row.names = FALSE, digits = 3)

## ===========================================================================
## Figures
## ===========================================================================
ggsave2 <- function(p, nm) ggsave(file.path(outdir, paste0("brms_compare_", nm, ".png")),
                                  p, width = 7, height = 4.5, dpi = 120)

# point estimates (study 1)
ggsave2(ggplot(res1$point, aes(merTools, brms)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.3, size = 0.8, color = "#1B9E77") + coord_equal() +
  labs(title = "Point estimates are essentially identical",
       subtitle = "Conditional mean per held-out student (STUDY 1)",
       x = "lme4 predict()", y = "brms posterior_epred()") +
  theme_minimal(base_size = 12), "point")

# coverage curves, both Gaussian studies x (seen/new)
cov_all <- do.call(rbind, list(
  transform(res1$cov_seen, study = "with meanses",    grp = "seen"),
  transform(res1$cov_new,  study = "with meanses",    grp = "new"),
  transform(res2$cov_seen, study = "without meanses", grp = "seen"),
  transform(res2$cov_new,  study = "without meanses", grp = "new")))
cov_long <- rbind(
  data.frame(cov_all[c("nominal", "study", "grp")], empirical = cov_all$merTools, method = "merTools"),
  data.frame(cov_all[c("nominal", "study", "grp")], empirical = cov_all$brms,     method = "brms"))
ggsave2(ggplot(cov_long, aes(nominal, empirical, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_line() + geom_point() + facet_grid(study ~ grp) +
  coord_equal(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
  labs(title = "Calibration: nominal vs. empirical coverage",
       x = "nominal level", y = "empirical coverage") +
  theme_minimal(base_size = 12) + theme(legend.position = "top"), "cover")

# GLMM probability agreement
ggsave2(ggplot(data.frame(merTools = mt_p, brms = brm_p), aes(merTools, brms)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.4, color = "#7570B3") + coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "GLMM: predicted P(tick) agrees with brms",
       subtitle = "grouseticks binomial model, held-out observations",
       x = "merTools predictInterval(type='probability')",
       y = "brms posterior_epred()") +
  theme_minimal(base_size = 12), "glmm")

cat("\nFigures written to ", outdir, "\n", sep = "")
