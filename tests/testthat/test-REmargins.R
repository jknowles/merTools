# Test REmargins
set.seed(51315)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

context("Test random effect marginalization works")

# mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:10,])
#
# test_that("Random effect marginals work for simple linear example", {
#   skip_on_travis()
#   skip_on_cran()
#   d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
#                    obs=1:100)
#   d$y <- simulate(~fac1+(1|grp),family = gaussian,
#                   newdata=d,
#                   newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
#                                  sigma = c(.23)), seed = 4548)[[1]]
#   subD <- d[sample(row.names(d), 1000),]
#
#   g1 <- lmer(y~fac1+(1|grp), data=subD)
#   d$fitted <- predict(g1, d)
#   #This suppresses the warning about no parallel backend registered
#   outs <- suppressWarnings(
#     predictInterval(g1, newdata = d, level = 0.9, n.sims = 1000,
#                     seed = 468,
#                     stat = 'mean', include.resid.var = TRUE)
#   )
#   outs <- cbind(d, outs); outs$coverage <- FALSE
#   outs$coverage <- outs$fitted <= outs$upr & outs$fitted >= outs$lwr
#   expect_true(all(outs$coverage))
#   expect_lt(abs(mean(outs$fit - outs$fitted)), .0005)
#   expect_lt(abs(mean(outs$fit - outs$y)), .01)
#   rm(outs)
# })
#
#
#
ggplot(out_w) + aes(x = obs, y = fit_Subject) +
  geom_line() +
  facet_wrap(~case)
