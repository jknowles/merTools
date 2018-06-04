#------------------------------------------------------------------------------
# ShinyMer
#-----------------------------------------------------------------------------

## Data and models
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
fm <- lmer( z ~ a*b + (1|r), data=d2)
fm2 <- lmer( z ~ a*b + (1+b|r), data=d2)


data(grouseticks)
grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
grouseticks$TICKS_BIN <- ifelse(grouseticks$TICKS >=1, 1, 0)
# GLMER 3 level + slope
form <- TICKS_BIN ~ YEAR + (1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
glmer3LevSlope  <- glmer(form, family="binomial",data=grouseticks,
                         control = glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun = 1e5)))

set.seed(3845)
d <- expand.grid(fac1=LETTERS[1:5], grp=factor(1:10), fac2 = LETTERS[10:20],
                 obs=1:25)
d$x <- runif(nrow(d))
d$y <- simulate(~ x + fac1 + fac2 + (1 + fac1|grp) + (1|obs), family = binomial,
                newdata=d,
                newparams=list(beta = rnorm(16),
                               theta = rnorm(16, 5, 1)))[[1]]
subD <- d[sample(row.names(d), 5000),]

g1 <- glmer(y ~ x + fac1 + fac2 + (1+fac1|grp) + (1|obs), data = subD, family = 'binomial')
g2 <- lmer(y ~ lectage + studage + (1+lectage|d) + (1|dept), data=InstEval)
# ----------------------------------------------------
# shinyMer blocks

shinyMer(fm)
shinyMer(fm2)
shinyMer(g1)
shinyMer(g2)
shinyMer(glmer3LevSlope)
shinyMer(glmer3LevSlope, simData = glmer3LevSlope@frame[22:25,])

