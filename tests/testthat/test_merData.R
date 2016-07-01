# -----------------------------------------------------
#-------------------------------------------------------

set.seed(51315)
library(lme4)
data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")


# Build out models
form <- TICKS ~ YEAR + HEIGHT +(1|BROOD) + (1|LOCATION) + (1|INDEX)
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
#Sanitize Names----
context("Sanitize Names")
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
#Strip Attributes----
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
#Random Observation----
context("Random observation")
################################################


test_that("A random observation can be sampled from a merMod", {
  data1 <- draw(glmer3Lev, type = 'random')
  data2 <- draw(lmerSlope2, type = 'random')
  data3 <- draw(lmerSlope1, type = 'random')
  data4 <- draw(glmer3LevSlope, type = 'random')
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

test_that("Random observation preserves factor levels", {
  data1 <- draw(glmer3Lev, type = 'random')
  data2 <- draw(lmerSlope2, type = 'random')
  data3 <- draw(lmerSlope1, type = 'random')
  data4 <- draw(glmer3LevSlope, type = 'random')
  expect_true(length(levels(data1$YEAR)) > length(unique(data1$YEAR)))
  expect_true(length(levels(data1$BROOD)) > length(unique(data1$BROOD)))
  expect_true(length(levels(data1$LOCATION)) > length(unique(data1$LOCATION)))
  expect_true(length(levels(data4$YEAR)) > length(unique(data4$YEAR)))
  expect_true(length(levels(data4$BROOD)) > length(unique(data4$BROOD)))
  expect_true(length(levels(data4$LOCATION)) > length(unique(data4$LOCATION)))
  # test levels are correct levels as well
})

###############################################
#Collapse frame----
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
#Super factor ----
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

test_that("SuperFactor handles new factor levels correctly", {
  fac1 <- factor(c("999", "888"))
  fac1a <- superFactor(fac1, fullLev = unique(grouseticks$BROOD))
  fac2 <- factor(c("Z16", "Z02", "Z05"))
  fac2a <- superFactor(fac2, fullLev = unique(Orthodont$Subject))
  expect_false(identical(levels(fac1), levels(fac1a)))
  expect_false(identical(levels(fac2), levels(fac2a)))
  expect_false(identical(levels(grouseticks$BROOD), levels(fac1a)))
  expect_false(identical(levels(Orthodont$Subject), levels(fac2a)))
  expect_equal(length(levels(fac1a)), length(levels(grouseticks$BROOD)) + 2)
  expect_equal(length(levels(fac2a)), length(levels(Orthodont$Subject)) + 3)
  expect_true(identical(levels(fac1a)[1:118], levels(grouseticks$BROOD)))
  expect_true(identical(levels(fac2a)[1:27], levels(Orthodont$Subject)))
})


###############################################
#Shuffle----
context("Shuffle")
################################################

test_that("Data can be shuffled", {
  expect_equal(nrow(Orthodont), nrow(merTools:::shuffle(Orthodont)))
  expect_equal(ncol(Orthodont), ncol(merTools:::shuffle(Orthodont)))
  expect_equal(nrow(grouseticks), nrow(merTools:::shuffle(grouseticks)))
  expect_equal(ncol(grouseticks), ncol(merTools:::shuffle(grouseticks)))
})

###############################################
#Find RE Quantiles----
context("Find RE Quantiles")
################################################

test_that("Errors and messages are issued", {
  expect_error(REquantile(glmer3Lev, 23, groupFctr = "BROOD"))
  expect_warning(REquantile(glmer3Lev, .23, groupFctr = "BROOD", term = "Cat"))
  expect_error(REquantile(glmer3Lev, .23, groupFctr = "Cat"))
  expect_error(REquantile(glmer3Lev, c(23, .56, .75), "BROOD"))
  expect_error(REquantile(glmer3Lev, c(.23, 56, .75), "BROOD"))
  expect_error(REquantile(glmer3Lev, c(.23, .56, 75), "BROOD"))
  expect_error(REquantile(glmer3Lev, c(.23, .56, 107), "BROOD"))
  expect_error(REquantile(glmer3Lev, c(-2, .56, .7), "BROOD"))
  expect_message(REquantile(lmerSlope1, .25, groupFctr = "Subject"))
  expect_warning(REquantile(lmerSlope2, c(.24), "Subject"))
  expect_warning(REquantile(lmerSlope2, c(.24), "Subject", term = "Cat"))
})

# what to do without intercepts (REquantile(lmerSlope2), c(.24), "Subject")
# test_that("Quantiles are returned correctly", {
#   myRE <- ranef(glmer3Lev)[["BROOD"]]
#   myRE <- myRE[order(myRE[, "(Intercept)"]), ,drop = FALSE]
#   rownames(myRE)[floor(23 / nrow(myRE)*100)]
#
#
# })

###############################################
#Test observation wiggle----
context("Test observation wiggle")
################################################

test_that("Row and column lengths are correct", {
  data1 <- grouseticks[5:9, ]
data1a <- wiggle(data1, var = "BROOD", values = c("606", "602", "537"))
  data1b <- wiggle(data1a, var = "YEAR", values = c("96", "97"))
  data2 <- grouseticks[3, ]
  data2a <- wiggle(data2, var = "BROOD", values = c("606", "602", "537"))
  data2b <- wiggle(data2a, var = "YEAR", values = c("96", "97"))
  data3 <- grouseticks[12:14, ]
  data3a <- wiggle(data3, var = "BROOD", values = c("606"))
  data3b <- wiggle(data3a, var = "YEAR", values = c("96", "97"))
  expect_equal(nrow(data1), 5)
  expect_equal(nrow(data1a), 15)
  expect_equal(nrow(data1b), 30)
  expect_equal(nrow(data2), 1)
  expect_equal(nrow(data2a), 3)
  expect_equal(nrow(data2b), 6)
  expect_equal(nrow(data3), 3)
  expect_equal(nrow(data3a), 3)
  expect_equal(nrow(data3b), 6)
  expect_equal(length(data1), length(data1a))
  expect_equal(length(data1a), length(data1b))
  expect_equal(length(data2), length(data2a))
  expect_equal(length(data2a), length(data2b))
  expect_equal(length(data3), length(data3a))
  expect_equal(length(data3a), length(data3b))
  data4 <- wiggle(data3, var = "BROOD",
                     values = REquantile(glmer3Lev,
                                         quantile = c(0.25, 0.5, 0.75),
                                         group = "BROOD"))
  expect_true(all(table(as.character(data4$BROOD),
                        as.character(data4$INDEX)) ==1))

})


test_that("Values are placed correctly", {
  data1 <- grouseticks[5:9, ]
  data1a <- wiggle(data1, var = "BROOD", values = c("606", "602", "537"))
  data1b <- wiggle(data1a, var = "YEAR", values = c("96", "97"))
  data2 <- grouseticks[3, ]
  data2a <- wiggle(data2, var = "BROOD", values = c("606", "602", "537"))
  data2b <- wiggle(data2a, var = "YEAR", values = c("96", "97"))
  data3 <- grouseticks[12:14, ]
  data3a <- wiggle(data3, var = "BROOD", values = c("606"))
  data3b <- wiggle(data3a, var = "YEAR", values = c("96", "97"))
  data4 <- Orthodont[15, ]
  data4a <- wiggle(data4, var = "age", values = c(10, 11, 12))
  data4b <- wiggle(data4a, var = "Sex", values = c("Male", "Female"))
  expect_false(any(unique(data1$BROOD) %in% unique(data1a$BROOD)))
  expect_false(any(unique(data1$BROOD) %in% unique(data1b$BROOD)))
  expect_false(any(unique(data1a$YEAR) %in% unique(data1b$YEAR)))
  expect_false(any(unique(data2$BROOD) %in% unique(data2a$BROOD)))
  expect_false(any(unique(data2$BROOD) %in% unique(data2b$BROOD)))
  expect_false(any(unique(data2a$YEAR) %in% unique(data2b$YEAR)))
  expect_false(any(unique(data3$BROOD) %in% unique(data3a$BROOD)))
  expect_false(any(unique(data3$BROOD) %in% unique(data3b$BROOD)))
  expect_false(any(unique(data3a$YEAR) %in% unique(data3b$YEAR)))
  expect_true(all(unique(data1a$BROOD) %in% c("606", "602", "537")))
  expect_true(all(unique(data1b$BROOD) %in% c("606", "602", "537")))
  expect_true(all(unique(data2a$BROOD) %in% c("606", "602", "537")))
  expect_true(all(unique(data2b$BROOD) %in% c("606", "602", "537")))
  expect_true(all(unique(data3a$BROOD) %in% c("606")))
  expect_true(all(unique(data3b$BROOD) %in% c("606")))
  expect_true(all(unique(data4a$age) %in% c(10, 11, 12)))
  expect_true(all(unique(data4b$age) %in% c(10, 11, 12)))
  expect_true(all(!unique(data1a$YEAR) %in% c("96", "97")))
  expect_true(all(unique(data1b$YEAR) %in% c("96", "97")))
  expect_true(all(!unique(data2a$YEAR) %in% c("96", "97")))
  expect_true(all(unique(data2b$YEAR) %in% c("96", "97")))
  expect_true(all(!unique(data3a$YEAR) %in% c("96", "97")))
  expect_true(all(unique(data3b$YEAR) %in% c("96", "97")))
  expect_true(all(unique(data4a$Sex) %in% c("Male", "Female")))
  expect_true(all(unique(data4b$Sex) %in% c("Male", "Female")))
})

###############################################
#Test average observation extraction----
context("Test average observation extraction")
################################################

test_that("Returns a single row", {
  data1 <- draw(glmer3Lev, type = 'average')
  data1a <- draw(glmer3LevSlope, type = 'average')
  data2 <- draw(lmerSlope1, type = 'average')
  expect_equal(nrow(data1), 1)
  expect_equal(nrow(data1a), 1)
  expect_equal(nrow(data2), 1)
})

test_that("Warnings and errors are correct", {
  expect_message(draw(lmerSlope1, type = 'average'))
  expect_warning(draw(lmerSlope2, type = 'average'))
  mylist2 <- list("YEAR" = "97", "LOCATION" = "16")
  expect_warning(draw(glmer3LevSlope, type = 'average', varList = mylist2))
  mylist3 <- list("YEAR" = "97", "LOCATION" = c("16", "56"))
  expect_warning(draw(glmer3LevSlope, type = 'average', varList = mylist3))
})

test_that("Subsets work", {
  mylist1 <- list("YEAR" = "97")
  data1 <- draw(glmer3LevSlope, type = 'average', varList = mylist1)
  data1a <- draw(glmer3LevSlope, type = 'average')
  expect_false(identical(data1, data1a))
  expect_equal(data1$TICKS, mean(grouseticks$TICKS[grouseticks$YEAR == "97"]))
  expect_equal(data1a$TICKS, mean(grouseticks$TICKS))
  mylist2 <- list("YEAR" = "97", "LOCATION" = "16")
  expect_warning(draw(glmer3LevSlope, type = 'average', varList = mylist2),
                 "less than 20 rows, averages may be problematic")
  mylist3 <- list("YEAR" = "97", "LOCATION" = c("16", "56"))
  expect_warning(draw(glmer3LevSlope, type = 'average', varList = mylist3),
                 "fewer than 3 rows, computing global average instead")
})

test_that("Nested specifications work", {
  library(ggplot2)
  mod1 <- lmer(sleep_total ~ bodywt + (1|vore/order), data=msleep)
  data1 <- draw(mod1, "random")
  expect_is(data1, "data.frame")
  data2 <- draw(mod1, "average")
  expect_is(data2, "data.frame")
  mylist1 <- list("vore" = "carni")
  mylist2 <- list("order" = "Cetacea")
  data1 <- draw(mod1, "random", varList = mylist1)
  expect_is(data1, "data.frame")
  expect_identical(as.character(data1$vore), "carni")
  data1 <- draw(mod1, "random", varList = mylist2)
  expect_is(data1, "data.frame")
  expect_identical(as.character(data1$order), "Cetacea")
  data1 <- suppressWarnings(draw(mod1, "average", varList = mylist1))
  expect_is(data1, "data.frame")
  expect_identical(as.character(data1$vore), "carni")
  data1 <- suppressWarnings(draw(mod1, "average", varList = mylist2))
  expect_is(data1, "data.frame")
  expect_identical(as.character(data1$order), "Cetacea")
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  data1 <- suppressWarnings(draw(fm1, type = "average", varList = list("Subject" = "308")))
  expect_is(data1, "data.frame")
  expect_identical(as.character(data1$Subject), "308")
})

