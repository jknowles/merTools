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

 #Build out model
 m1.form <- Reaction ~ Days + (Days | Subject)
 m1.df <- sleepstudy
 m1.new.df <- sleepstudy

 #Estimate model
 m1  <- lmer(m1.form, data=m1.df)
 display(m1)

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
   data.frame(
     median = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
     lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
     upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))),
     colNAs = apply(merBoot$t, 2, function(x) sum(is.na(x)))
   )
 }

#Step 2: Our method ####
 source("F:/RSandbox/merTools/R/merExtract.R")

 predictInterval <- function(model, newdata, level = 0.95){
   require(mvtnorm)
   require(lme4)
   reTerms <- names(ngrps(model))
   if(length(reTerms) > 1){
     stop("Multiple grouping terms not yet implemented.")
   }
   group <- reTerms[1]
   outs <- data.frame(fit = rep(NA, nrow(newdata)),
                      lwr = rep(NA, nrow(newdata)),
                      upr = rep(NA, nrow(newdata)))

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

   total <- nrow(newdata)
   pb <- txtProgressBar(min=0, max=total, style=3)

   #Begin cleaning up the loop
   betaSim <- rmvnorm(1000, mean = fixef(model), sigma = as.matrix(vcov(model)))
   newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X

   for(i in 1:nrow(newdata)){
     caseID <- newdata[i, group]
     reSim <- rmvnorm(1000,
                      mean = reMeansExtract(model, group = group, case=caseID),
                      sigma = reMatrixExtract(model, group = group, case = caseID))
     totCoef <- combineREFE(reSim, betaSim)
     #this gets the model matrix our of the merMod ... this will only work if newdata == originaldata
     yhat <- totCoef %*% newdata.modelMatrix[i,]
     outs[i, "fit"] <- mean(yhat)
     outs[i, "upr"] <- as.numeric(quantile(yhat, 1 - ((1-level)/2)))
     outs[i, "lwr"] <- as.numeric(quantile(yhat, ((1-level)/2)))
     setTxtProgressBar(pb,i)
   }
   close(pb)
   return(outs)
 }

 kf.time <- system.time(
   kf.method <- predictInterval(m1, m1.new.df)
 )

#Step 3: Summarize and Compare Results####
 ##Summarize and Combine
 boot1.sum <- sumBoot(boot1.m1)
 boot1.sum$id <- 1:nrow(boot1.sum)
 boot1.sum <- select(boot1.sum, id, fit=median, lwr, upr) %>%
   gather(stat, Boot.1, fit, lwr, upr)

 boot2.sum <- sumBoot(boot2.m1)
 boot2.sum$id <- 1:nrow(boot2.sum)
 boot2.sum <- select(boot2.sum, id, fit=median, lwr, upr) %>%
   gather(stat, Boot.2, fit, lwr, upr)

 boot3.sum <- sumBoot(boot3.m1)
 boot3.sum$id <- 1:nrow(boot3.sum)
 boot3.sum <- select(boot3.sum, id, fit=median, lwr, upr) %>%
   gather(stat, Boot.3, fit, lwr, upr)

 kf.method$id <- 1:nrow(kf.method)
 kf.method <- select(kf.method, id, fit, lwr, upr) %>%
   gather(stat, KF, fit, lwr, upr)

 eval <- merge(boot1.sum, boot2.sum, by=c("id", "stat"))
 eval <- merge(eval, boot3.sum, by=c("id","stat"))
 eval <- merge(eval, kf.method, by=c("id","stat"))

 eval <- eval %>%
   gather(type, value, Boot.1, Boot.2, Boot.3, KF) %>%
   spread(stat, value)

 eval$id[eval$type=="Boot.1"] <- eval$id[eval$type=="Boot.1"] - .15
 eval$id[eval$type=="Boot.2"] <- eval$id[eval$type=="Boot.2"] - .05
 eval$id[eval$type=="Boot.3"] <- eval$id[eval$type=="Boot.3"] + .05
 eval$id[eval$type=="KF"]     <- eval$id[eval$type=="KF"]     + .15

 #Gen both point predictions from predict.merMod for reference
 m1.new.df$id <- 1:nrow(m1.new.df)
 m1.new.df$with.u <- predict(m1, newdata=m1.new.df, re.form=NULL)
 m1.new.df$no.u <- predict(m1, newdata=m1.new.df, re.form=NA)
 m1.new.df <- gather(m1.new.df, ptype, fit, with.u, no.u) %>%
   arrange(id)

 #Diagnostic plots
 ggplot(aes(y=fit,x=id, color=type), data=eval[1:180,]) +
   geom_point(size=I(3)) +
   geom_linerange(aes(ymax=upr, ymin=lwr), size=I(1)) +
   geom_point(aes(shape=ptype), data=m1.new.df[1:90,], color="black") +
   theme_bw() +
   labs(x="Index", y="Prediction Interval", title="95% Prediction interval by method for select individuals")

 ggplot(aes(x=fit), data=eval) +
   geom_density(aes(color=type)) +
   geom_vline(x=mean(test$y)) +
   geom_density(data=summarize(group_by(test, d), fit = mean(y)), color="black") +
   theme_bw() +
   labs(x = "Estimate", y="Density",
        title="Average point estimates by prediction type \n (vertical black line is sample grand mean, \n black density is distribution of sample group means)")




