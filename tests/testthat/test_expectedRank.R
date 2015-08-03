# Test expected rank

#Using 2 of sample models from test_merExtract.R

set.seed(51315)
library(lme4)

# Sleepstudy
m1  <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)

m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Wackier example
data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
m3 <- lmer(distance ~ age + (0 + age + nsex|Subject), data=Orthodont)

# two grouping factors
data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")

form <- TICKS ~ YEAR + HEIGHT +(HEIGHT|BROOD) + (1|INDEX)
m4  <- glmer(form, family="poisson",data=grouseticks,
                     control = glmerControl(optimizer="Nelder_Mead",
                                            optCtrl=list(maxfun = 1e5)))

form <- TICKS ~ YEAR + HEIGHT + (0 + HEIGHT|BROOD) + (1|INDEX)
m5  <- glmer(form, family="poisson",data=grouseticks,
                     control = glmerControl(optimizer="Nelder_Mead",
                                            optCtrl=list(maxfun = 1e5)))

#Custom Expectation Functions
  expect_correct_dim <- function(merMod, factor=NULL, term=NULL) {
    if (is.null(factor))
      n.levels <- nrow(ranef(merMod)[[1]])
    else
      n.levels <- nrow(ranef(merMod)[[factor]])
    ER <- expectedRank(merMod, factor, term)
    expect_true(nrow(ER)==n.levels &&
                ncol(ER)==5 &&
                colnames(ER)[2:5]==c("theta", "varTheta", "ER", "pctER") &&
                class(ER)=="data.frame")
  }


###############################################
context("Testing expected rank")
###############################################
test_that("expectedRank parameters work and dont work as intended", {
  expect_correct_dim(m1)
  expect_correct_dim(m1, factor="Subject")
  expect_correct_dim(m1, term="(Intercept)")
  expect_correct_dim(m1, factor="Subject", term="(Intercept)")

  expect_correct_dim(m2, term="(Intercept)")
  expect_correct_dim(m2, term="Days")

  expect_correct_dim(m3, factor="Subject", term="age")
  expect_correct_dim(m3, factor="Subject", term="nsex")

  expect_correct_dim(m4, factor="BROOD", term="(Intercept)")
  expect_correct_dim(m4, factor="INDEX", term="(Intercept)")

  expect_correct_dim(m5, factor="BROOD")
  expect_correct_dim(m5, factor="INDEX")

  expect_error(expectedRank(m4), "Must specify which grouping factor when there are more than one")
  expect_error(expectedRank(m4, factor="BROOD"), "Must specify which random coefficient when there are more than one per selected grouping factor")
  expect_error(expectedRank(m3, factor="Subject"), "Must specify which random coefficient when there are more than one per selected grouping factor")
  expect_error(expectedRank(m3, term="int"), "undefined columns selected")
})

test_that("Ranks have the correct range", {
  numGrps <- nrow(ranef(m1)[[1]])
  expect_true(max(expectedRank(m1)$ER) <= numGrps)
  expect_true(min(expectedRank(m1)$ER) >= 1)
  expect_equal(cor(expectedRank(m1)$ER, rank(ranef(m1)[[1]])), -0.09582297)
})

test_that("Percentile ranks have the correct range", {
   expect_true(max(expectedRank(m1)$pctER) <= 100)
   expect_true(min(expectedRank(m1)$pctER) >= 0)
})



