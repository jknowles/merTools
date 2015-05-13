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

 #Cannonical Models for testing purposes
 ##1) Sleepstudy eg - lmer with random slope and random intercept
   m1.form <- Reaction ~ Days + (Days | Subject)
   m1.df <- sleepstudy
   m1.new.df <- m1.df
 ##2) Verbal Aggresion eg - lmer (could to glmer(logit)) with 2 levels
   m2.form <- r2 ~ (Anger + Gender + btype + situ)^2 + (1|id) + (1|item)
   m2.df <- VerbAgg
   m2.new.df <- m2.df

#Step 1: Our method ####
 predictInterval <- function(model, newdata, level = 0.95){
   #Depends
   require(mvtnorm)
   require(lme4)
   #Prep Output
   outs <- newdata

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
   outs$fit <- apply(yhat,1,function(x) as.numeric(quantile(x, .5)))
   outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2))))
   outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2))))
   #Close it out
   return(outs)
 }

#Step 2: Unit Test Function####
  predictInterval.test <- function(model.form, model.df, predict.df, idvar=NULL, ...) {
   require(lme4); require(dplyr); require(tidyr); require(ggplot2);
   ##Estimate model
   modelEstimation.time <- system.time(
     m1  <- lmer(model.form, data=model.df, ...)
   )
   ##If it does not have one, add unique identifier to predict.df
   if (is.null(idvar)) {
     predict.df$.newID <- paste("newID", rownames(predict.df), sep="")
   }
   ##Functions for bootMer() and objects
   ####Return predicted values from bootstrap
   mySumm <- function(.) {
     predict(., newdata=predict.df, re.form=NULL, type="link")
   }
   ####Collapse bootstrap into median, 95% PI
   sumBoot <- function(merBoot) {
     data.frame(merBoot$data,
                median = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
                lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
                upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))),
                colNAs = apply(merBoot$t, 2, function(x) sum(is.na(x)))
     )
   }
   ##Bootstrap
   ####Method 1: parametric, re-estimate BLUPS
   cat("\n Bootstrap Method 1")
   boot1.time <- system.time(
     boot1 <- bootMer(m1, mySumm, nsim=1000,
                         use.u=FALSE, type="parametric",
                         .progress = "txt", PBargs=list(style=3))
   )
   ####Method 2: parametric, conditional on estimated BLUPS
   cat("\n Bootstrap Method 2")
   boot2.time <- system.time(
     boot2 <- bootMer(m1, mySumm, nsim=1000,
                         use.u=TRUE, type="parametric",
                         .progress = "txt", PBargs=list(style=3))
   )
   ####Method 3: semiparametric (draw from resid), conditional on estimated BLUPS
   cat("\n Bootstrap Method 3")
   boot3.time <- system.time(
     boot3 <- bootMer(m1, mySumm, nsim=1000,
                         use.u=TRUE, type="semiparametric",
                         .progress = "txt", PBargs=list(style=3))
   )
   ##Our Method
   kf.time <- system.time(
     kf.method <- predictInterval(m1, predict.df)
   )
   ##Compare Times
   compare.time <- rbind(modelEstimation.time, kf.time, boot1.time, boot2.time, boot3.time)

   ##Summarize and compare results
   boot1.sum <- sumBoot(boot1)
   boot2.sum <- sumBoot(boot2)
   boot3.sum <- sumBoot(boot3)

   eval <- merge(predict.df, boot1.sum) %>%
     rename(boot1.fit=median, boot1.lwr=lwr, boot1.upr=upr)
   eval <- merge(eval, boot2.sum) %>%
     rename(boot2.fit=median, boot2.lwr=lwr, boot2.upr=upr)
   eval <- merge(eval, boot3.sum) %>%
     rename(boot3.fit=median, boot3.lwr=lwr, boot3.upr=upr)
   eval <- merge(eval, kf.method) %>%
     rename(KF.fit=fit, KF.lwr=lwr, KF.upr=upr)

   #Check if nrow(eval) still equals nrow(predict.df) because it should
   if (nrow(eval)!=nrow(predict.df)) {
     stop("Something happened when merging bootstrap summaries together ...")
   }

   ##Add lmer yhats (predict.merMod) on there
   eval$with.u <- predict(m1, newdata=predict.df, re.form=NULL, type="link")
   eval$no.u <- predict(m1, newdata=predict.df, re.form=NA, type="link")

   ##Diagnostic plots
   ####Data prep
   Eval <- sample_n(eval, 30) %>%
     gather(var, value, starts_with("boot"), starts_with("KF")) %>%
     mutate(
       simMethod = sub("[[:punct:]][0-9A-Za-z]*", "", var),
       stat = sub("[0-9A-Za-z]*[[:punct:]]", "", var)
     ) %>%
     select(-var) %>%
     spread(stat, value) %>%
     group_by(.newID) %>%
     arrange(simMethod) %>%
     mutate(
       x=row_number(simMethod),
       x=(x-mean(x))/10
     ) %>%
     group_by(simMethod) %>%
     mutate(
       x=row_number(.newID)+x
     )


   ####Direct Comparison Plot
   p1 <- ggplot(aes(y=fit,x=x, color=simMethod), data=Eval) +
           geom_point(size=I(3)) +
           geom_linerange(aes(ymax=upr, ymin=lwr), size=I(1)) +
           geom_point(shape="with.u", color="black", size=I(4),
                      data=summarize(group_by(Eval, .newID), x=mean(x), fit=mean(with.u))) +
           geom_point(shape="no.u", color="black", size=I(4),
                      data=summarize(group_by(Eval, .newID), x=mean(x), fit=mean(no.u))) +
           theme_bw() +
           labs(x="Index", y="Prediction Interval", title="95% Prediction interval by method for 30 random obs")

   #####Distribution of fitted values
   p2 <- ggplot(aes(x=fit), data=Eval) +
           geom_density(aes(color=Sim.Method)) +
           geom_vline(x=mean(m1@resp$y)) +
           geom_density(data=summarize(
                               group_by_(
                                 cbind(y=m1@resp$y, as.data.frame(m1@flist)),
                                 names(m1@flist)),
                              fit = mean(y)), color="black") +
           theme_bw() +
           labs(x = "Estimate", y="Density",
                title="Average point estimates by prediction type \n (vertical black line is sample grand mean, \n black density is distribution of sample group means)")


   ##Close it out
   return(
     list(
       compareTimes=compare.time,
       bootstraps = list(boot1, boot2, boot3),
       kf.method = kf.method,
       model=m1,
       evalData = eval,
       plot.directCompare = p1,
       plot.pointEstimateDist = p2
       )
     )
 }

#Step 3: Summarize and Compare Results####
 #debug(predictInterval.test)
 cannonical.1 <-  predictInterval.test(m1.form, m1.df, m1.new.df)
 cannonical.1$plot.directCompare

 cannonical.2 <-  predictInterval.test(m2.form, m2.df, m2.new.df)
 cannonical.2$plot.directCompare

 save.image()



