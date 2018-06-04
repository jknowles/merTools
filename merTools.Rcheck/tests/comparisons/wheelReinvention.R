#What exactly are the values spit out by fitted(sim(merMod), merMod)???
#Does it produce what we are trying to produce with predictInterval?

#The answer is no because it doesn't re-sample the RFX (I am fairly certain)
#The following does show the arm::sim() side-by-side with our method.  It looks
#like the effect of the differences is that our within-group regression lines are
#closer to the population average regression line.

library(arm); library(abind); library(mvtnorm);
set.seed(95371)

data(sleepstudy)

model <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)

#Take 5 simulations
test.sim <- sim(model, 1000)

#Display simulated coefs
#coef(test.sim)
ranef.array <- coef(test.sim)$ranef[[1]]
fixef.matrix <- coef(test.sim)$fixef

#Expand and permute fixefs to conform to ranefs for elementwise addition
dim(ranef.array)
dim(fixef.matrix)

#For now do it by hand
fixef.array <- array(data = rep(fixef.matrix, dim(ranef.array)[2]),
                     dim  = c(dim(fixef.matrix), dim(ranef.array)[2]))
fixef.array <- aperm(fixef.array, c(1,3,2))
dim(fixef.array)

combo.array <- fixef.array+ranef.array

#Extract model.matrix for our own multiplication
model.matrix <- lFormula(model@call, data=model@frame)$X
dim(model.matrix)
expanded.combo.array <- combo.array[,model@flist[[1]],]
expanded.combo.array <- aperm(expanded.combo.array, c(3,2,1))

#Matrix multiplication with arrays ???

##The following yeilds a 180x180 matrix for each prediction because I have
##multiplied ALL obs by ALL possible group coefficients
myCalc.large <- abind(
                  lapply(1:dim(expanded.combo.array)[3],
                         function(i) model.matrix %*% expanded.combo.array[,,i]),
                  along=3)

##I only need the diagonal of each 180x180 matrix so we get the value of
##the values for obs i and the coefficients for obs i
myCalc <- abind(lapply(1:dim(myCalc.large)[3],
                       function(x) diag(myCalc.large[,,x])), along=2)

isTRUE(all.equal(myCalc,fitted(test.sim, model), check.attributes=FALSE))

##Compare with predictInterval()
checkPI <- predictInterval(model, sleepstudy, nsim = 1000, predict="link")$yhat

myCalc.lwr <- apply(myCalc,1,function(x) as.numeric(quantile(x, .025)))
myCalc.fit <- apply(myCalc,1,function(x) as.numeric(quantile(x, .500)))
myCalc.upr <- apply(myCalc,1,function(x) as.numeric(quantile(x, .975)))

checkPI.lwr <- apply(checkPI,1,function(x) as.numeric(quantile(x, .025)))
checkPI.fit <- apply(checkPI,1,function(x) as.numeric(quantile(x, .500)))
checkPI.upr <- apply(checkPI,1,function(x) as.numeric(quantile(x, .975)))

plot.data <- rbind(
               data.frame(model="Arm.sim",
                          x=(1:180)-.15,
                          lwr=myCalc.lwr,
                          fit=myCalc.fit,
                          upr=myCalc.upr),
               data.frame(model="predictInterval",
                          x=(1:180)+.15,
                          lwr=checkPI.lwr,
                          fit=checkPI.fit,
                          upr=checkPI.upr))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=model, position="dodge"), data=plot.data) +
geom_point() +
geom_linerange()

library(dplyr)
calc.sigma <- sqrt(1/rgamma(1000, 0.5*lme4:::df.residual.merMod(model), 0.5*getME(model, "devcomp")$cmp[["pwrss"]]))
sim.sigma <- test.sim@sigma
data.frame(
  model=c(rep("ARM",1000),rep("Carl",1000)),
  sigma=c(sim.sigma, calc.sigma)
  ) %>%
  qplot(sigma, color=model, data=., geom="density")
