# Test helper functions
# Trimming data frame----
context("Trimming data frame")
test_that("Trimming results in correct size", {
  data(InstEval)
  trimDat <- merTools:::trimModelFrame(InstEval)
  expect_gt(nrow(InstEval), nrow( merTools:::trimModelFrame(InstEval)))
  expect_equal(nrow(trimDat), 4065)
  cbpp$obs <- 1:nrow(cbpp)
  d1 <- cbpp
  d1$y <- d1$incidence / d1$size
  gm2 <- glmer(y ~ period +
                  (1 | herd),
                family = binomial, data = d1, nAGQ = 9, weights = d1$size)
  trimDat <- merTools:::trimModelFrame(gm2@frame)
  expect_is(trimDat, "data.frame")
  expect_equal(nrow(trimDat), 18)
})

test_that("Trimming does not corrupt order", {
  tmp <- InstEval[1:10, ]
  trimDat <- merTools:::trimModelFrame(InstEval)
  trimDat <- rbind(tmp, trimDat)
  expect_lt(nrow(trimDat), nrow(tmp) + nrow(InstEval))
  row.names(tmp) <- NULL
  row.names(trimDat) <- NULL
  expect_identical(tmp, trimDat[1:10, ])
})

# subBoot and Theta----
context("subBoot and Theta")

test_that("Can extract theta from a fit model", {
  set.seed(404)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)
  d$y <- simulate(~fac1+(1|grp),family = gaussian,
                  newdata=d,
                  newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                                 sigma = c(.23)))[[1]]
  subD <- d[sample(row.names(d), 1000),]
  g1 <- lmer(y~fac1+(1|grp), data=subD)

  g1b <- lm(y ~ fac1, data = subD)

  expect_equal(thetaExtract(g1), 0.2085, tolerance = .05)
  expect_error(thetaExtract(g1b))

  z1 <- subBoot(g1, 500, FUN = thetaExtract, R = 10)
  expect_is(z1, "data.frame")
  expect_equal(nrow(z1), 11)
  expect_equal(ncol(z1), 2)
})

# Test formula Build-----
context("Test formula build")

test_that("Formula works for additive functions", {
  n <- 20
  x <- y <- rnorm(n)
  z <- rnorm(n)
  r <- sample(1:5, size=n, replace=TRUE)
  d <- data.frame(x,y,z,r)
  d2 <- expand.grid(a=factor(1:4),b=factor(1:4),rep=1:10)
  n <- nrow(d2)
  d2 <- transform(d2,r=sample(1:5, size=n, replace=TRUE),
                  z=rnorm(n))
  d2 <- subset(d2,!(a=="4" & b=="4"))
  fm <- lmer( z ~ a + b + (1|r), data=d2)
  expect_is(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a + b"))
})


test_that("Formula works for interactions", {
  n <- 200
  x <- y <- rnorm(n)
  z <- rnorm(n)
  r <- sample(1:5, size=n, replace=TRUE)
  d <- data.frame(x,y,z,r)
  d2 <- expand.grid(a=factor(1:4),b=factor(1:4), c = factor(1:4), rep=1:10)
  n <- nrow(d2)
  d2 <- transform(d2,r=sample(1:5, size=n, replace=TRUE),
                  z=rnorm(n))
  d2 <- subset(d2,!(a=="4" & b=="4"))
  d2$x <- rnorm(nrow(d2))
  fm <- lmer( z ~ a * b + c +  (1|r), data=d2)
  expect_is(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b + c"))
  fm <- lmer( z ~ a * b * c +  (1|r), data=d2)
  expect_is(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b * c"))
  fm <- lmer( z ~ a * b * c + x + I(x^2) + (1 + c|r), data=d2)
  expect_is(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b * c + x + I(x^2)"))
})

