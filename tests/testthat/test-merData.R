# -----------------------------------------------------
# Test framework includes tests for multiple intercepts and
# multiple slopes to ensure extraction of random effects
# works in these scenarios
#-------------------------------------------------------

set.seed(51315)
library(lme4)
data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")


# Build out models
form <- TICKS ~ YEAR + HEIGHT +(1|BROOD) + (1|INDEX) + (1|LOCATION)
glmer3Lev  <- glmer(form, family="poisson",data=grouseticks,
                    control = glmerControl(optimizer="Nelder_Mead",
                                           optCtrl=list(maxfun = 1e5)))
# GLMER 3 level + slope
form <- TICKS ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
glmer3LevSlope  <- glmer(form, family="poisson",data=grouseticks,
                         control = glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun = 1e5)))

# GLMER 2 level
# data(VerbAgg)
# fmVA <- glmer(r2 ~ Anger + Gender + btype + situ +
#                 (1|id) + (1|item), family = binomial, data =
#                 VerbAgg)

# Sleepstudy
lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Wackier example
data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
lmerSlope2 <- lmer(distance ~ age + (0 + age + nsex|Subject), data=Orthodont)



###############################################
context("Santize Names")
################################################

test_that("Sanitize names renames variables in data.frame", {
  badMod <- lmer(distance ~ factor(Sex) + (0 + age + nsex|Subject),
                 data=Orthodont)
  expect_false(identical(names(badMod@frame),
                         names(merTools:::sanitizeNames(badMod@frame))))
  expect_is(merTools:::sanitizeNames(badMod@frame), "data.frame")
  expect_identical(names(merTools:::sanitizeNames(badMod@frame))[2], "Sex")
  expect_identical(names(badMod@frame)[2], "factor(Sex)")
})


###############################################
context("Strip attributes")
################################################

test_that("Attributes can be stripped from data.frame", {
  full <- names(attributes(lmerSlope1@frame))
  redu <- names(attributes(merTools:::stripAttributes(lmerSlope1@frame)))
  redu2 <- names(attributes(merTools:::stripAttributes(glmer3LevSlope@frame)))
  expect_true(length(full) > length(redu))
  expect_true(all(redu %in% full))
  expect_true(all(redu %in% c("names", "row.names", "class")))
  expect_true(all(redu2 %in% c("names", "row.names", "class")))
})


###############################################
context("Random observation")
################################################

randomObs <- function(merMod){
  out <- merMod@frame[sample(1:nrow(merMod@frame), 1),]
  chars <- !sapply(out, is.numeric)
  for(i in names(out[, chars])){
    out[, i] <- superFactor(out[, i], fullLev = unique(merMod@frame[, i]))
  }
  out <- stripAttributes(out)
  return(out)
}

test_that("A random observation can be sampled from a merMod", {
  data1 <- randomObs(glmer3Lev)
  data2 <- randomObs(lmerSlope2)
  data3 <- randomObs(lmerSlope1)
  data4 <- randomObs(glmer3LevSlope)
  expect_equal(nrow(data1), 1)
  expect_equal(nrow(data2), 1)
  expect_equal(nrow(data3), 1)
  expect_equal(nrow(data4), 1)
  expect_equal(ncol(data1), 6)
  expect_equal(ncol(data2), 4)
  expect_equal(ncol(data3), 3)
  expect_equal(ncol(data4), 6)
  expect_identical(names(data1), names(glmer3Lev@frame))
  expect_identical(names(data2), names(lmerSlope2@frame))
  expect_identical(names(data3), names(lmerSlope1@frame))
  expect_identical(names(data4), names(glmer3LevSlope@frame))
  expect_false(identical(names(attributes(data1)), names(attributes(glmer3Lev@frame))))
  expect_false(identical(names(attributes(data2)), names(attributes(lmerSlope2@frame))))
  expect_false(identical(names(attributes(data3)), names(attributes(lmerSlope1@frame))))
  expect_false(identical(names(attributes(data4)), names(attributes(glmer3LevSlope@frame))))
  expect_false("formula" %in% names(attributes(data1)))
  expect_false("formula" %in% names(attributes(data2)))
  expect_false("formula" %in% names(attributes(data3)))
  expect_false("formula" %in% names(attributes(data4)))
})

###############################################
context("Collapse frame")
################################################

test_that("Collapsing a dataframe results in single row", {
  data1 <- merTools:::collapseFrame(Orthodont)
  data2 <- merTools:::collapseFrame(grouseticks)
  expect_equal(length(data1), length(Orthodont))
  expect_equal(length(data2), length(grouseticks))
  expect_equal(nrow(data1), 1)
  expect_equal(nrow(data2), 1)
  expect_equal(data1$distance, mean(Orthodont$distance))
  expect_equal(data1$distance, mean(Orthodont$distance))
  expect_equal(data1$age, mean(Orthodont$age))
  expect_equal(data1$nsex, mean(Orthodont$nsex))
  expect_equal(data1$nsexage, mean(Orthodont$nsexage))
  expect_equal(data2$TICKS, mean(grouseticks$TICKS))
  expect_equal(data2$HEIGHT, mean(grouseticks$HEIGHT))
  expect_equal(data2$cHEIGHT, mean(grouseticks$cHEIGHT))
  expect_equal(data2$meanTICKS, mean(grouseticks$meanTICKS))
})

###############################################
context("Subset by a list")
################################################

test_that("Data can be subset by a list", {
  list11 <- list("Sex" = "Male")
  list12 <- list("Sex" = "Male", "Subject" = "M05")
  data11 <- merTools:::subsetList(Orthodont, list11)
  data12 <- merTools:::subsetList(Orthodont, list12)
  list21 <- list("YEAR" = "95")
  list22 <- list("LOCATION" = "32", "BROOD" = "503")
  data21 <- merTools:::subsetList(grouseticks, list21)
  data22 <- merTools:::subsetList(grouseticks, list22)
  expect_equal(length(data11), length(Orthodont))
  expect_equal(length(data21), length(grouseticks))
  expect_equal(length(data12), length(Orthodont))
  expect_equal(length(data22), length(grouseticks))
  expect_equal(nrow(data11), 64)
  expect_equal(nrow(data21), 117)
  expect_equal(nrow(data12), 4)
  expect_equal(nrow(data22), 0)
})

###############################################
context("Super factor")
################################################

test_that("Unobserved factor levels can be respected", {
  fac1 <- factor(c("502", "503"))
  fac1a <- superFactor(fac1, fullLev = unique(grouseticks$BROOD))
  fac2 <- factor(c("M16", "M02", "M05"))
  fac2a <- superFactor(fac2, fullLev = unique(Orthodont$Subject))
  expect_false(identical(levels(fac1), levels(fac1a)))
  expect_false(identical(levels(fac2), levels(fac2a)))
  expect_true(identical(levels(grouseticks$BROOD), levels(fac1a)))
  expect_true(identical(levels(Orthodont$Subject), levels(fac2a)))
  expect_equal(length(levels(fac1a)), 118)
  expect_equal(length(levels(fac2a)), 27)
})

###############################################
context("Shuffle")
################################################

test_that("Data can be shuffled", {
  expect_equal(nrow(Orthodont), nrow(merTools:::shuffle(Orthodont)))
  expect_equal(ncol(Orthodont), ncol(merTools:::shuffle(Orthodont)))
  expect_equal(nrow(grouseticks), nrow(merTools:::shuffle(grouseticks)))
  expect_equal(ncol(grouseticks), ncol(merTools:::shuffle(grouseticks)))
})

