# test merModList functions

#Do merModList objects get built and work----
context("Do merModList objects get built and work")

old_warn <- getOption("warn")
options(warn = -1)

test_that("simple cases work", {
  # skip_on_cran()
  library(blme)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  split <- sample(x = LETTERS[9:15], size = nrow(d), replace=TRUE)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  out <- split(d, split)
  rm(split)
  # TODO change tolerances
  g1 <- lmerModList(formula = y~fac1+(1|grp), data=out,
                    control= lmerControl(check.conv.grad = .makeCC("warning", tol= 2e-3)))
  expect_is(g1, "merModList")
  g2 <- blmerModList(formula = y~fac1+(1|grp), data=out,
                     control= lmerControl(check.conv.grad = .makeCC("warning", tol= 2e-3)))
  expect_is(g2, "merModList")
  expect_false(class(g1[[1]]) == class(g2[[1]]))

  split <- sample(x = LETTERS[1:20], size = nrow(InstEval), replace=TRUE)
  out <- split(InstEval, split)
  rm(split)
  g1 <- lmerModList(formula = y ~ lectage + studage + (1|d) + (1|dept),
                    data=out,
                    control= lmerControl(check.conv.grad = .makeCC("warning", tol = 2e-4)))
  expect_is(g1, "merModList")

})

test_that("print methods work for merModList", {
  # skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  split <- sample(x = LETTERS[9:15], size = nrow(d), replace=TRUE)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  out <- split(d, split)
  rm(split);
  g1 <- lmerModList(formula = y~fac1+(1|grp), data=out,
                    control= lmerControl(check.conv.grad = .makeCC("warning", tol= 1e-2)));
  {sink("NUL"); zz <- print(g1);
    sink()}
  expect_is(zz, "list")
  zz <- summary(g1)
  expect_is(zz, "summary.merModList")

})

# Numerical accuracy of merModList print method----
context("Numerical accuracy of merModList print method")

test_that("print method for merModList works in general case", {
  # skip_on_cran()
  data(grouseticks)
  grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
  grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
  grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
  form <- TICKS_BIN ~ HEIGHT +(1 + HEIGHT|BROOD) + (1|YEAR)
  modDat <- vector(5, mode="list")
  for(i in 1:length(modDat)){
    modDat[[i]] <- grouseticks[sample(x=1:nrow(grouseticks), size = nrow(grouseticks),
                                      replace=FALSE),]
  }

  g1 <- glmerModList(formula = form,
                    data = modDat, family = "binomial",
                    control = glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun = 1e6),
                                           check.conv.grad = .makeCC("warning", tol= 1e-2)))
  g1T <- glmer(form, family = "binomial", data = grouseticks,
               control = glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun = 1e6),
                                      check.conv.grad = .makeCC("warning", tol= 1e-2)))

  expect_equal(VarCorr(g1)$stddev$BROOD, attr(VarCorr(g1T)$BROOD, "stddev"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$stddev$YEAR, attr(VarCorr(g1T)$YEAR, "stddev"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$correlation$BROOD, attr(VarCorr(g1T)$BROOD, "corre"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$correlation$YEAR, attr(VarCorr(g1T)$YEAR, "corre"),
               tolerance = 0.0001)

  form <- TICKS_BIN ~ HEIGHT +(1|BROOD)
  g1 <- glmerModList(formula = form,
                     data = modDat, family = "binomial",
                     control = glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun = 1e6),
                                            check.conv.grad = .makeCC("warning", tol= 1e-2)))
  g1T <- glmer(form, family = "binomial", data = grouseticks,
               control = glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun = 1e6),
                                      check.conv.grad = .makeCC("warning", tol= 1e-2)))

  expect_equal(VarCorr(g1)$stddev$BROOD, attr(VarCorr(g1T)$BROOD, "stddev"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$stddev$YEAR, attr(VarCorr(g1T)$YEAR, "stddev"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$correlation$BROOD, attr(VarCorr(g1T)$BROOD, "corre"),
               tolerance = 0.0001)
  expect_equal(VarCorr(g1)$correlation$YEAR, attr(VarCorr(g1T)$YEAR, "corre"),
               tolerance = 0.0001)

})

#ICC function----
context("ICC function")

test_that("ICC function works", {
  ICC1 <- ICC(outcome = "Reaction", group = "Subject", data = sleepstudy)
  expect_is(ICC1, "numeric")
  expect_equal(ICC1, 0.3948896, tol = .001)
})

options(warn= old_warn)
