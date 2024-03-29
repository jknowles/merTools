# Test helper functions
set.seed(51315)
# Trimming data frame----

test_that("Trimming results in correct size", {
  skip_on_cran()
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
  expect_s3_class(trimDat, "data.frame")
  expect_equal(nrow(trimDat), 18)
})

test_that("Trimming does not corrupt order", {
  skip_on_cran()
  tmp <- InstEval[1:10, ]
  trimDat <- merTools:::trimModelFrame(InstEval)
  trimDat <- rbind(tmp, trimDat)
  expect_lt(nrow(trimDat), nrow(tmp) + nrow(InstEval))
  row.names(tmp) <- NULL
  row.names(trimDat) <- NULL
  expect_identical(tmp, trimDat[1:10, ])
})

# subBoot and Theta----
# context("subBoot and Theta")

test_that("Can extract theta from a fit model", {
  skip_on_cran()
  set.seed(404)
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:100)

  suppressMessages({
  d$y <- simulate( ~ fac1 + (1|grp),
                  newdata = d,
                  family = gaussian,
                  newparams = list(
                    "theta" = 0.22,
                    "beta" = c(2,1,3,4,7),
                    "sigma" = 0.23))[[1]]
  })

  subD <- d[sample(row.names(d), 1000),]
  g1 <- lmer(y~fac1+(1|grp), data=subD)

  g1b <- lm(y ~ fac1, data = subD)

  expect_equal(thetaExtract(g1), 0.2285, tolerance = 0.1)
  expect_error(thetaExtract(g1b))

  z1 <- suppressMessages({
    subBoot(g1, 500, FUN = thetaExtract, R = 10)
  })
  expect_s3_class(z1, "data.frame")
  expect_equal(nrow(z1), 11)
  expect_equal(ncol(z1), 2)
})

# Test formula Build-----
# context("Test formula build")

test_that("Formula works for additive functions", {
  skip_on_cran()
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
  expect_s3_class(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a + b"))
})


test_that("Formula works for interactions", {
  skip_on_cran()
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
  suppressMessages({
    fm <- lmer( z ~ a * b + c +  (1|r), data=d2)
  })

  expect_s3_class(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b + c"))
  suppressMessages({
    fm <- lmer( z ~ a * b * c +  (1|r), data=d2)
  })

  expect_s3_class(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b * c"))
  suppressMessages({
    fm <- lmer( z ~ a * b * c + x + I(x^2) + (1 + c|r), data=d2)
  })

  expect_s3_class(merTools:::formulaBuild(fm), "formula")
  expect_identical(merTools:::formulaBuild(fm), as.formula("z ~ a * b * c + x + I(x^2)"))
})


test_that("Build model matrix produces matrices of the right size", {
  skip_on_cran()
  d <- expand.grid(fac1 = LETTERS[1:5],
                   grp = letters[11:20],
                   obs = 1:50)
  suppressMessages({
    d$y <- simulate(~fac1 + (1 | grp), family = binomial,
                    newdata = d,
                    newparams = list( "theta" = c(.33),
                                      "beta" = c(2,-1,3,-2,1.2)),
                    seed =634)[[1]]

  })

  subD <- d[sample(row.names(d), 1200), ]

  g1 <- glmer(y~fac1+(1|grp), data=subD, family = 'binomial')
  d$fitted <- predict(g1, d)

  mm <- merTools:::buildModelMatrix(g1, newdata = d, which = "full")
  expect_true(inherits(mm, "matrix") || inherits(mm, "Matrix"))
  expect_equal(dim(mm), c(2500, 15))

}
)
