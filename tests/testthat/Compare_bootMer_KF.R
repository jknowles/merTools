#This file compares the prediction interval from the "Correct" bootmer method
#to our quick and dirty method to see how they differ.

#I renamed some things to ease my typing

#Prep R####
 #rm(list=ls())
 library(lme4)
 library(arm)
 library(mvtnorm)
 library(dplyr)
 library(tidyr)
 library(ggplot2)

 set.seed(51315)
 data(InstEval)
 data(sleepstudy)
 data(VerbAgg)

 #Build out model
 m1.form <- Reaction ~ Days + (Days | Subject)
 m1.df <- sleepstudy
 m1.new.df <- m1.df

 #Estimate model
 m1.estimation.time <- system.time(
   m1  <- lmer(m1.form, data=m1.df)
 )
 print(m1)

#Step 1: Get bootMer... gold standard ####
 ##Define function for parameter extraction (inspired from bootMer help file)
 ##Pulling:
 ##  Individual level predicted values conditional on ALL random effects
 mySumm <- function(.) {
    predict(., newdata=m1.new.df, re.form=NULL, type="link")
    }

 ##Maybe try all three (valid) ways listed in bootMer doc in the order listed in details
 boot1.time <- system.time(
   boot1.m1 <- bootMer(m1, mySumm, nsim=1000,
                        use.u=FALSE, type="parametric",
                        .progress = "txt", PBargs=list(style=3))
 )

 boot2.time <- system.time(
   boot2.m1 <- bootMer(m1, mySumm, nsim=1000,
                        use.u=TRUE, type="parametric",
                        .progress = "txt", PBargs=list(style=3))
 )

 boot3.time <- system.time(
   boot3.m1 <- bootMer(m1, mySumm, nsim=1000,
                        use.u=TRUE, type="semiparametric",
                        .progress = "txt", PBargs=list(style=3))
 )

 sumBoot <- function(merBoot) {
   data.frame(merBoot$data,
     median = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
     lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
     upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))),
     colNAs = apply(merBoot$t, 2, function(x) sum(is.na(x)))
   )
 }

#Step 2: Our method ####
 ##"sub" functions
 source("F:/RSandbox/merTools/R/merExtract.R")

 reMatrixExtract <- function(model, group, case){
   slotNum <- ranef(model)[[group]]
   slotNum$seq <- 1:nrow(slotNum)
   slotNum <- slotNum[case, ]$seq
   mat <- as.matrix(attr(ranef(model, condVar=TRUE)[[group]],
                         which = "postVar")[, , slotNum])
   row.names(mat) <- names(ranef(model)[[group]])
   colnames(mat) <- row.names(mat)
   return(mat)
 }

 reMeansExtract <- function(model, group, case){
   rerow <- ranef(model)[[group]][case, ]
   out <- as.numeric(rerow)
   names(out) <- names(ranef(model)[[group]])
   return(out)
 }

 combineREFE <- function(reSim, betaSim){
   comboCols <- intersect(colnames(reSim), colnames(betaSim))
   REonly <- setdiff(colnames(reSim), colnames(betaSim))
   FEonly <- setdiff(colnames(betaSim), colnames(reSim))
   totCoef <- cbind(reSim[, comboCols] + betaSim[, comboCols])
   colnames(totCoef) <- comboCols
   totCoef <- cbind(totCoef, reSim[, REonly])
   totCoef <- cbind(totCoef, betaSim[, FEonly])
   return(totCoef)
 }

 ##Main Function
 predictInterval <- function(model, newdata, level = 0.95){
   #Depends
   require(mvtnorm)
   require(lme4)
   #Prep Output
   outs <- data.frame(fit = rep(NA, nrow(newdata)),
                      lwr = rep(NA, nrow(newdata)),
                      upr = rep(NA, nrow(newdata)))

   #Sort out all the levels
   reTerms <- names(ngrps(model))
   n.reTerms = length(reTerms)

   ##The following 3 lines produce a matrix of linear predictors created from the fixed coefs
   betaSim <- rmvnorm(1000, mean = fixef(model), sigma = as.matrix(vcov(model)))
   newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X
   fixed.xb <- newdata.modelMatrix %*% t(betaSim)

   ##Draw from random effects distributions for each level and merge onto newdata
   reSim <- NULL
   nSims = 1000
   for (j in seq_along(reTerms)) {
     group=reTerms[j]
     reMeans <- array(ranef(model)[[group]])
     reMatrix <- attr(ranef(model, condVar=TRUE)[[group]], which = "postVar")
     reSim[[group]] <- data.frame(rownames(reMeans), matrix(NA, nrow=nrow(reMeans), ncol=nSims))
     colnames(reSim[[group]]) <- c(group, paste("sim", 1:nSims, sep=""))
     for (k in 1:nrow(reMeans)) {
       lvl = rownames(reMeans)[k]
       reSim[[group]][k,2:ncol(reSim[[group]])] <- rmvnorm(nSims, mean=as.matrix(reMeans[k,]), sigma=as.matrix(reMatrix[,,k]))
     }
     cnames <- colnames(reSim[[group]])
     reSim[[group]] <- merge(newdata, reSim[[group]], by=group, all.x=TRUE)
     reSim[[group]] <- as.matrix(reSim[[group]][,setdiff(cnames, group)])
   }

   #Calculate yhat as sum of components
   yhat <- fixed.xb + apply(simplify2array(reSim), c(1,2), sum)

   #Output prediction intervals
   outs$fit <- apply(yhat,1,mean)
   outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2))))
   outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2))))
   #Close it out
   return(outs)
 }

 kf.time <- system.time(
   kf.method <- predictInterval(m1, m1.new.df)
 )

 ##Compare Times
 compare.time <- rbind(kf.time, boot1.time, boot2.time, boot3.time)
 rm(kf.time, boot1.time, boot2.time, boot3.time)

#Step 3: Summarize and Compare Results####
 ##Summarize and Combine
 boot1.sum <- sumBoot(boot1.m1)
 boot1.sum$.new.id <- 1:nrow(boot1.sum)
 boot1.sum <- select(boot1.sum, .new.id, fit=median, lwr, upr) %>%
   gather(stat, Boot.1, fit, lwr, upr)

 boot2.sum <- sumBoot(boot2.m1)
 boot2.sum$.new.id <- 1:nrow(boot2.sum)
 boot2.sum <- select(boot2.sum, .new.id, fit=median, lwr, upr) %>%
   gather(stat, Boot.2, fit, lwr, upr)

 boot3.sum <- sumBoot(boot3.m1)
 boot3.sum$.new.id <- 1:nrow(boot3.sum)
 boot3.sum <- select(boot3.sum, .new.id, fit=median, lwr, upr) %>%
   gather(stat, Boot.3, fit, lwr, upr)

 kf.method$.new.id <- 1:nrow(kf.method)
 kf.method <- select(kf.method, .new.id, fit, lwr, upr) %>%
   gather(stat, KF, fit, lwr, upr)

 eval <- merge(boot1.sum, boot2.sum, by=c(".new.id", "stat"))
 eval <- merge(eval, boot3.sum, by=c(".new.id","stat"))
 eval <- merge(eval, kf.method, by=c(".new.id","stat"))

 eval <- eval %>%
   gather(Sim.Method, value, Boot.1, Boot.2, Boot.3, KF) %>%
   spread(stat, value)

 eval$.new.id[eval$Sim.Method=="Boot.1"] <- eval$.new.id[eval$Sim.Method=="Boot.1"] - .15
 eval$.new.id[eval$Sim.Method=="Boot.2"] <- eval$.new.id[eval$Sim.Method=="Boot.2"] - .05
 eval$.new.id[eval$Sim.Method=="Boot.3"] <- eval$.new.id[eval$Sim.Method=="Boot.3"] + .05
 eval$.new.id[eval$Sim.Method=="KF"]     <- eval$.new.id[eval$Sim.Method=="KF"]     + .15

 #Gen both point predictions from predict.merMod for reference
 m1.new.df$.new.id <- 1:nrow(m1.new.df)
 m1.new.df$with.u <- predict(m1, newdata=m1.new.df, re.form=NULL, type="link")
 m1.new.df$no.u <- predict(m1, newdata=m1.new.df, re.form=NA, type="link")
 m1.new.df <- gather(m1.new.df, pred.merMod.type, fit, with.u, no.u) %>%
   arrange(.new.id)

 #Diagnostic plots
 ggplot(aes(y=fit,x=.new.id, color=Sim.Method), data=eval[1:100,]) +
   geom_point(size=I(3)) +
   geom_linerange(aes(ymax=upr, ymin=lwr), size=I(1)) +
   geom_point(aes(shape=pred.merMod.type), data=m1.new.df[1:50,], color="black") +
   theme_bw() +
   labs(x="Index", y="Prediction Interval", title="95% Prediction interval by method for select individuals")


 ggplot(aes(x=fit), data=eval) +
   geom_density(aes(color=type)) +
   geom_vline(x=mean(test$y)) +
   geom_density(data=summarize(group_by(test, d), fit = mean(y)), color="black") +
   theme_bw() +
   labs(x = "Estimate", y="Density",
        title="Average point estimates by prediction type \n (vertical black line is sample grand mean, \n black density is distribution of sample group means)")




