library(microbenchmark)

# ClassFilter <- function(x) inherits(get(x), 'lm' ) & !inherits(get(x), 'gl
set.seed(101)


# Small

lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

bench1 <- microbenchmark(
  predict(lmerSlope1, newdata = sleepstudy[1:100,]),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 100, stat = 'median',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 100, stat = 'mean',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 100, stat = 'median',
                  include.resid.var = FALSE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 1000, stat = 'median',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.8, nsim = 1000, stat = 'median',
                  include.resid.var = TRUE),
  times = 10
)

bench2 <- microbenchmark(
  predict(lmerSlope1, newdata = sleepstudy[1:100,]),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 100, stat = 'median',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 200, stat = 'mean',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 400, stat = 'median',
                  include.resid.var = FALSE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 800, stat = 'median',
                  include.resid.var = TRUE),
  predictInterval(lmerSlope1, newdata = sleepstudy[1:100,],
                  level = 0.9, nsim = 1600, stat = 'median',
                  include.resid.var = TRUE),
  times = 10
)

# Medium
d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
                 obs=1:100)
d$y <- simulate(~fac1+(1|grp),family = gaussian,
                newdata=d,
                newparams=list(beta=c(2,1,3,4,7), theta=c(.25),
                               sigma = c(.23)))[[1]]
g1 <- lmer(y~fac1+(1|grp), data=d)


bench3 <- microbenchmark(predictInterval(g1, newdata = d[1:100, ], level = 0.9,
                               nsim = 50,
                               stat = 'mean', include.resid.var = TRUE),
               predictInterval(g1, newdata = d[1:200, ], level = 0.9,
                               nsim = 50,
                               stat = 'mean', include.resid.var = TRUE),
               predictInterval(g1, newdata = d[1:400, ], level = 0.9,
                               nsim = 50,
                               stat = 'mean', include.resid.var = TRUE),
               predictInterval(g1, newdata = d[1:800, ], level = 0.9,
                               nsim = 50,
                               stat = 'mean', include.resid.var = TRUE),
               times = 10)


# Large

g2 <- lmer(y ~ lectage + studage + (1+lectage|d) + (1|dept), data=InstEval)
d2 <- InstEval[1:1000, ]
outs1a <- predictInterval(g2, newdata = d2, level = 0.8, n.sims = 500,
                          stat = 'mean', include.resid.var=TRUE)


bench4 <- microbenchmark(predictInterval(g2, newdata = d2[1:100, ], level = 0.9,
                                         nsim = 500,
                                         stat = 'mean', include.resid.var = TRUE),
                         predictInterval(g2, newdata = d2[1:200, ], level = 0.9,
                                         nsim = 500,
                                         stat = 'mean', include.resid.var = TRUE),
                         predictInterval(g2, newdata = d2[1:400, ], level = 0.9,
                                         nsim = 500,
                                         stat = 'mean', include.resid.var = TRUE),
                         predictInterval(g2, newdata = d2[1:800, ], level = 0.9,
                                         nsim = 500,
                                         stat = 'mean', include.resid.var = TRUE),
                         times = 10)

# set.seed(101)
# d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10),
#                  obs=1:50)
# d$y <- simulate(~fac1+(1|grp),family = binomial,
#                 newdata=d,
#                 newparams=list(beta=c(2,-1,3,-2,1.2), theta=c(.33)))[[1]]
# subD <- d[sample(row.names(d), 1200),]
# g1 <- glmer(y~fac1+(1|grp), data=subD, family = 'binomial')
# d$fitted <- predict(g1, d)
#
#
# outs <- predictInterval(g1, newdata = d, level = 0.95, nsim = 500,
#                         stat = 'mean', include.resid.var = FALSE,
#                         type = 'linear.prediction')
#
#
# g2 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
# d1 <- InstEval[1:100, ]
#






