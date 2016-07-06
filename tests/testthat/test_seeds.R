# 'seed' options in draw, REsim, FEsim, predictInterval and subBoot----
context("'seed' options in draw, REsim, FEsim, predictInterval and subBoot")

test_that("Equivalent seeds return equivalent results", {
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

  d1a <- draw(fm1, type="random", seed=1234)
  d2  <- draw(fm1, type="random", seed=456)
  d1b <- draw(fm1, type="random", seed=1234)

  r1a <- REsim(fm1, 25, seed=1234)
  r2  <- REsim(fm1, 25, seed=456)
  r1b <- REsim(fm1, 25, seed=1234)

  f1a <- FEsim(fm1, 25, seed=1234)
  f2  <- FEsim(fm1, 25, seed=456)
  f1b <- FEsim(fm1, 25, seed=1234)

  p1a <- predictInterval(fm1, newdata=sleepstudy[1:10,], seed=1234)
  p2  <- predictInterval(fm1, newdata=sleepstudy[1:10,], seed=456)
  p1b <- predictInterval(fm1, newdata=sleepstudy[1:10,], seed=1234)

  s1a <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20, seed=1234)
  s2  <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20, seed=456)
  s1b <- subBoot(fm1, n = 160, FUN = thetaExtract, R = 20, seed=1234)

  expect_identical(d1a, d1b)
  expect_identical(r1a, r1b)
  expect_identical(f1a, f1b)
  expect_identical(p1a, p1b)
  expect_identical(s1a, s1b)

  expect_false(identical(d1a, d2))
  expect_false(identical(r1a, r2))
  expect_false(identical(f1a, f2))
  expect_false(identical(p1a, p2))
  expect_false(identical(s1a, s2))
})
