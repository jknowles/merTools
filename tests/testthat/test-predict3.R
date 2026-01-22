# Special cases - rank deficiency----

test_that("Prediction intervals are accurate with interaction terms and rank deficiency", {
  skip_on_ci()
  skip_on_cran()
  set.seed(54656)
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

  fm <- lmer( z ~ a*b + (1|r), data=d2) |> suppressMessages()
  expect_s3_class(predictInterval(fm, newdata = d2[1:10, ]), "data.frame")

  newPred <- predictInterval(fm, newdata = d2, level = 0.8, n.sims = 1500,
                             stat = 'median', include.resid.var = FALSE,
                             seed = 2342)
  truPred <- predict(fm, newdata = d2)
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/15)
  fm2 <- lmer( z ~ a*b + (1+b|r), data=d2) |> suppressMessages()
  #In the call below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  newPred <- suppressWarnings(predictInterval(fm2, newdata = d2, level = 0.8,
                                              n.sims = 1000,
                                              stat = 'median', include.resid.var = FALSE))
  truPred <- predict(fm2, newdata = d2)
  expect_s3_class(newPred, "data.frame")
  expect_equal(mean(newPred$fit - truPred), 0, tolerance = sd(truPred)/10)
})

# Test the simResults----


test_that("simResults option behaves", {
  skip_on_cran()
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  preds1 <- predictInterval(m1, newdata = sleepstudy[1:5, ])
  preds2 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE)
  expect_null(attr(preds1, "sim.results"))
  expect_true(is.matrix(attr(preds2, "sim.results")))
  out <- attr(preds2, "sim.results")
  expect_equal(ncol(out), 1000)
  expect_equal(nrow(out), 5)
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  preds1 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE,
                            which = "random", seed = 23151)
  preds2 <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "fixed",
                            returnSims = TRUE, seed = 23151)
  preds3 <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "all",
                            returnSims = TRUE, seed = 23151)
  preds4 <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                            returnSims = TRUE, seed = 23151)
  preds1b <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                             returnSims = TRUE,
                             which = "random", seed = 23151,
                             include.resid.var = FALSE)
  preds2b <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "fixed",
                             returnSims = TRUE, seed = 23151,
                             include.resid.var = FALSE)
  preds3b <- predictInterval(m1, newdata = sleepstudy[1:5, ], which = "all",
                             returnSims = TRUE, seed = 23151,
                             include.resid.var = FALSE)
  preds4b <- predictInterval(m1, newdata = sleepstudy[1:5, ],
                             returnSims = TRUE, seed = 23151,
                             include.resid.var = FALSE)
  expect_true(is.matrix(attr(preds1, "sim.results")))
  expect_gt(abs(mean(attr(preds1, "sim.results") - attr(preds2, "sim.results"))),
            200)
  expect_true(is.matrix(attr(preds2, "sim.results")))
  expect_gt(abs(mean(attr(preds4, "sim.results") - attr(preds1, "sim.results"))),
            200)
  expect_gt(abs(mean(attr(preds4, "sim.results") - attr(preds1, "sim.results"))),
            abs(mean(attr(preds1, "sim.results") - attr(preds2, "sim.results"))))
  expect_type(attr(preds3, "sim.results"), "list")
  expect_true(is.list(attr(preds3, "sim.results")))
  expect_true(is.matrix(attr(preds1b, "sim.results")))
  expect_true(is.matrix(attr(preds2b, "sim.results")))
  expect_true(is.list(attr(preds3b, "sim.results")))
  expect_true(is.matrix(attr(preds4b, "sim.results")))
  # Check that samples are wider for include.resid.var = TRUE
  expect_gt(quantile(attr(preds1, "sim.results"), probs = 0.9) - quantile(attr(preds1b, "sim.results"), probs = 0.9),
            20)
  expect_lt(quantile(attr(preds1, "sim.results"), probs = 0.1) - quantile(attr(preds1b, "sim.results"), probs = 0.1),
            -20)
  expect_gt(quantile(attr(preds2, "sim.results"), probs = 0.9) - quantile(attr(preds2b, "sim.results"), probs = 0.9),
            20)
  expect_lt(quantile(attr(preds2, "sim.results"), probs = 0.1) - quantile(attr(preds2b, "sim.results"), probs = 0.1),
            -20)
  expect_gt(quantile(attr(preds4, "sim.results"), probs = 0.9) - quantile(attr(preds4b, "sim.results"), probs = 0.9),
            15)
  expect_lt(quantile(attr(preds4, "sim.results"), probs = 0.1) - quantile(attr(preds4b, "sim.results"), probs = 0.1),
            -15)
})

# Test out of sample predictions----

test_that("predictInterval makes predictions without observed outcome", {
  skip_on_ci()
  skip_on_cran()
  possNames <- expand.grid(letters,LETTERS)
  possNames <- paste(possNames[, 1], possNames[, 2])
  newFac <- sample(possNames, 32)
  modData <- data.frame(
    y = rnorm(500),
    x = rnorm(500),
    team_name = sample(newFac, 500, replace = TRUE)
  )
  modData$y[251:500] <- rep(NA, 250)
  m0 <- lmer(y ~ x + (1|team_name), data = modData[1:250,]) |> suppressMessages()
  #In the calls below we are getting warnings because our call to mvtnorm::rmvnorm
  #is shotting a warning when mean and sigma of multivariate distribution are
  #zero using the method="chol
  testPreds1 <- suppressWarnings(predictInterval(m0, newdata = modData[, c(3, 2, 1)]))
  testPreds2 <- suppressWarnings(predictInterval(m0, newdata = modData[1:250, c(2, 3, 1)]))
  testPreds3 <- suppressWarnings(predictInterval(m0, newdata = modData[251:500,]))
  expect_s3_class(testPreds1, "data.frame")
  expect_s3_class(testPreds2, "data.frame")
  expect_s3_class(testPreds3, "data.frame")
})

# Input validation checks----



test_that("dplyr objects are successfully coerced", {
  skip_on_cran()
  set.seed(101)
  data(sleepstudy)
  m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
  predData <- sleepstudy |>
    dplyr::group_by(Subject) |>
    dplyr::summarise(Days = mean(Days))
  expect_warning(predictInterval(m1, newdata = predData),
                 regexp = "newdata is tbl_df or tbl object from dplyr package")
  #Suppress the warning that we tested for above
  preds2 <- suppressWarnings(predictInterval(m1, newdata = predData, n.sims=2000))
  expect_s3_class(preds2, "data.frame")
  predData2 <- as.data.frame(predData)
  preds1 <- predictInterval(m1, newdata = predData2, n.sims=2000)
  expect_true(sum(preds1$fit - preds2$fit) > -50 & sum(preds1$fit - preds2$fit) < 50)
})

# Model type warnings for non-binomial GLMM----


test_that("Warnings issued", {
  skip_on_cran()
  d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                   obs=1:50)
  d$y <- simulate(~fac1+(1|grp),family = poisson,
                  newdata=d,
                  newparams=list(beta=c(2,-1,3,-2,1.2), theta=c(.33)),
                  seed = 5636)[[1]] |> suppressMessages()
  g1 <- glmer(y~fac1+(1|grp), data=d, family = 'poisson')
  expect_warning(predictInterval(g1, newdata = d[1:100,]),
                 regexp = "Prediction for NLMMs or GLMMs that are not mixed")
  rm(list = ls())
})
