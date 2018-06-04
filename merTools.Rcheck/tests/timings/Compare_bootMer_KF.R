# #This file compares the prediction interval from the "Correct" bootmer method
# #to our quick and dirty method to see how they differ.
#
# #I renamed some things to ease my typing
#
# #Prep R####
#  #rm(list=ls())
#  library(lme4)
#  library(arm)
#  library(mvtnorm)
#  library(dplyr)
#  library(tidyr)
#  library(ggplot2)
#  library(knitr)
#  library(RPushbullet)
#
#  set.seed(51315)
#  data(InstEval)
#  data(sleepstudy)
#  data(VerbAgg)
#
#  #Cannonical Models for testing purposes
#  ##1) Sleepstudy eg - lmer with random slope and random intercept
#    m1.form <- Reaction ~ Days + (Days | Subject)
#    m1.df <- sleepstudy
#    m1.new.df <- m1.df
#  ##2) Verbal Aggression eg - lmer (could to glmer(logit)) with 2 levels
#    m2.form <- r2 ~ (Anger + Gender + btype + situ)^2 + (1|id) + (1|item)
#    m2.df <- VerbAgg
#    m2.new.df <- m2.df
#  ##3) Verbal Aggression eg - glmer(logit) with 2 levels
#    m3.form <- r2 ~ (Anger + Gender + btype + situ)^2 + (1|id) + (1|item)
#    m3.df <- VerbAgg
#    m3.new.df <- m3.df
#  ##4) Instructer Evaluation et - lmer cross-classified, with two different random slopes
#    m4.form <- y ~ studage + lectage + service + dept + (studage|s) + (lectage|d)
#    m4.df <- InstEval
#    m4.new.df <- m4.df
#
# #Step 1: Our method ####
# #  predictInterval <- function(model, newdata, level = 0.95, nsim=1000, stat="median", predict.type="link"){
# #    #Depends
# #    require(mvtnorm)
# #    require(lme4)
# #    #Prep Output
# #    outs <- newdata
# #
# #    #Sort out all the levels
# #    reTerms <- names(ngrps(model))
# #    n.reTerms = length(reTerms)
# #
# #    ##The following 3 lines produce a matrix of linear predictors created from the fixed coefs
# #    betaSim <- rmvnorm(nsim, mean = fixef(model), sigma = as.matrix(vcov(model)))
# #    newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X
# #    fixed.xb <- newdata.modelMatrix %*% t(betaSim)
# #
# #    ##Draw from random effects distributions for each level and merge onto newdata
# #    reSim <- NULL
# #    for (j in seq_along(reTerms)) {
# #      group=reTerms[j]
# #      reMeans <- array(ranef(model)[[group]])
# #      reMatrix <- attr(ranef(model, condVar=TRUE)[[group]], which = "postVar")
# #      reSim[[group]] <- data.frame(rownames(reMeans), matrix(NA, nrow=nrow(reMeans), ncol=nsim))
# #      colnames(reSim[[group]]) <- c(group, paste("sim", 1:nsim, sep=""))
# #      for (k in 1:nrow(reMeans)) {
# #        lvl = rownames(reMeans)[k]
# #        reSim[[group]][k,2:ncol(reSim[[group]])] <- rmvnorm(nsim, mean=as.matrix(reMeans[k,]), sigma=as.matrix(reMatrix[,,k]))
# #      }
# #      cnames <- colnames(reSim[[group]])
# #      reSim[[group]] <- merge(newdata, reSim[[group]], by=group, all.x=TRUE)
# #      reSim[[group]] <- as.matrix(reSim[[group]][,setdiff(cnames, group)])
# #    }
# #
# #    #Calculate yhat as sum of components
# #    yhat <- fixed.xb + apply(simplify2array(reSim), c(1,2), sum)
# #
# #    #Output prediction intervals
# #    if (stat=="median") {
# #      outs$fit <- apply(yhat,1,function(x) as.numeric(quantile(x, .5)))
# #    }
# #    if (stat=="mean") {
# #      outs$fit <- apply(yhat,1,mean)
# #    }
# #    outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2))))
# #    outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2))))
# #    if (predict.type=="response") {
# #      outs$fit <- model@resp$family$linkinv(outs$fit)
# #      outs$upr <- model@resp$family$linkinv(outs$upr)
# #      outr$lwr <- model@resp$family$linkinc(outs$lwr)
# #    }
# #    #Close it out
# #    return(outs)
# #  }
#
# #Step 2: Unit Test Function####
#  predictInterval.test <- function(model.form, model.df, model.type="lmer",
#                                   predict.df, predict.type="link", stat="median",
#                                   idvar=NULL, nSims=1000, ...) {
#    require(lme4); require(dplyr); require(tidyr); require(ggplot2);
#    ##Estimate model
#    if (model.type=="lmer") {
#      modelEstimation.time <- system.time(
#        m1  <- lmer(model.form, data=model.df, ...)
#      )
#    }
#    if (model.type=="glmer") {
#      modelEstimation.time <- system.time(
#        m1  <- glmer(model.form, data=model.df, ...)
#      )
#    }
#    if (model.type=="blmer") {
#      modelEstimation.time <- system.time(
#        m1  <- blmer(model.form, data=model.df, ...)
#      )
#    }
#    if (model.type=="bglmer") {
#      modelEstimation.time <- system.time(
#        m1  <- bglmer(model.form, data=model.df, ...)
#      )
#    }
#    ##If it does not have one, add unique identifier to predict.df
#    if (is.null(idvar)) {
#      predict.df$.newID <- paste("newID", rownames(predict.df), sep="")
#    }
#    ##Functions for bootMer() and objects
#    ####Return predicted values from bootstrap
#    mySumm <- function(.) {
#      predict(., newdata=predict.df, re.form=NULL, type=predict.type)
#    }
#    ####Collapse bootstrap into median, 95% PI
#    sumBoot <- function(merBoot, stat) {
#      if (stat=="median") {
#         fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))
#      }
#      if (stat=="mean") {
#         fit = apply(merBoot$t, 2, function(x) mean(x, na.rm=TRUE))
#      }
#      return(
#        data.frame(merBoot$data,
#                   fit = fit,
#                   lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
#                   upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))),
#                   colNAs = apply(merBoot$t, 2, function(x) sum(is.na(x)))
#        )
#      )
#    }
#    ##Bootstrap
#    ####Method 1: parametric, re-estimate BLUPS
#    cat("\n Bootstrap Method 1")
#    boot1.time <- system.time(
#      boot1 <- bootMer(m1, mySumm, nsim=nSims,
#                          use.u=FALSE, type="parametric",
#                          .progress = "txt", PBargs=list(style=3))
#    )
#    ####Method 2: parametric, conditional on estimated BLUPS
#    cat("\n Bootstrap Method 2")
#    boot2.time <- system.time(
#      boot2 <- bootMer(m1, mySumm, nsim=nSims,
#                          use.u=TRUE, type="parametric",
#                          .progress = "txt", PBargs=list(style=3))
#    )
#    ####Method 3: semiparametric (draw from resid), conditional on estimated BLUPS
#    cat("\n Bootstrap Method 3")
#    boot3.time <- system.time(
#      boot3 <- bootMer(m1, mySumm, nsim=nSims,
#                          use.u=TRUE, type="semiparametric",
#                          .progress = "txt", PBargs=list(style=3))
#    )
#    ##Our Method
#    kf.time <- system.time(
#      kf.method <- predictInterval(m1, predict.df)
#    )
#    ##Compare Times
#    compare.time <- rbind(modelEstimation.time, kf.time, boot1.time, boot2.time, boot3.time)
#
#    ##Summarize and compare results
#    boot1.sum <- sumBoot(boot1, stat=stat)
#    boot2.sum <- sumBoot(boot2, stat=stat)
#    boot3.sum <- sumBoot(boot3, stat=stat)
#
#    eval <- merge(predict.df, boot1.sum) %>%
#      rename(boot1.fit=fit, boot1.lwr=lwr, boot1.upr=upr)
#    eval <- merge(eval, boot2.sum) %>%
#      rename(boot2.fit=fit, boot2.lwr=lwr, boot2.upr=upr)
#    eval <- merge(eval, boot3.sum) %>%
#      rename(boot3.fit=fit, boot3.lwr=lwr, boot3.upr=upr)
#    eval <- merge(eval, kf.method) %>%
#      rename(KF.fit=fit, KF.lwr=lwr, KF.upr=upr)
#
#    #Check if nrow(eval) still equals nrow(predict.df) because it should
#    if (nrow(eval)!=nrow(predict.df)) {
#      stop("Something happened when merging bootstrap summaries together ...")
#    }
#
#    ##Add lmer yhats (predict.merMod) on there
#    eval$with.u <- predict(m1, newdata=predict.df, re.form=NULL, type=predict.type)
#    eval$no.u <- predict(m1, newdata=predict.df, re.form=NA, type=predict.type)
#
#    ##Create summary statistics
#    piCoveragePct <- function(ref.upr, ref.lwr, new.upr, new.lwr) {
#      pct <- ifelse(ref.upr < new.lwr | ref.lwr > new.upr, 0,
#             ifelse(ref.upr < new.upr & ref.lwr > new.lwr, 1,
#             ifelse(ref.upr > new.upr & ref.lwr < new.lwr, (new.upr-new.lwr)/(ref.upr-ref.lwr),
#             ifelse(ref.upr < new.upr, (ref.upr-new.lwr)/(ref.upr-ref.lwr),
#             ifelse(ref.lwr > new.lwr, (new.upr-ref.lwr)/(ref.upr-ref.lwr), NA)))))
#      return(pct)
#    }
#    ###Pct of other PI covered by KF PI
#    eval$boot1.coverage <- piCoveragePct(eval$boot1.upr, eval$boot1.lwr, eval$KF.upr, eval$KF.lwr)
#    eval$boot2.coverage <- piCoveragePct(eval$boot2.upr, eval$boot2.lwr, eval$KF.upr, eval$KF.lwr)
#    eval$boot3.coverage <- piCoveragePct(eval$boot3.upr, eval$boot3.lwr, eval$KF.upr, eval$KF.lwr)
#    ###Does KF PI contain point estimates
#    eval$contain.boot1  <- piCoveragePct(eval$boot1.fit, eval$boot1.fit, eval$KF.upr, eval$KF.lwr)
#    eval$contain.boot2  <- piCoveragePct(eval$boot2.fit, eval$boot2.fit, eval$KF.upr, eval$KF.lwr)
#    eval$contain.boot3  <- piCoveragePct(eval$boot3.fit, eval$boot3.fit, eval$KF.upr, eval$KF.lwr)
#    eval$contain.with.u <- piCoveragePct(eval$with.u   , eval$with.u   , eval$KF.upr, eval$KF.lwr)
#    eval$contain.no.u   <- piCoveragePct(eval$no.u     , eval$no.u     , eval$KF.upr, eval$KF.lwr)
#    ###Point estimate "bias"
#    eval$distance.boot1  <- eval$KF.fit - eval$boot1.fit
#    eval$distance.boot2  <- eval$KF.fit - eval$boot2.fit
#    eval$distance.boot3  <- eval$KF.fit - eval$boot3.fit
#    eval$distance.with.u <- eval$KF.fit - eval$with.u
#    eval$distance.no.u   <- eval$KF.fit - eval$no.u
#
#    ##Close it out
#    return(
#      list(
#        compareTimes=compare.time,
#        bootstraps = list(boot1, boot2, boot3),
#        kf.method = kf.method,
#        model=m1,
#        evalData = eval
#        )
#      )
#  }
#
# #Step 2b: Post processing predictInterval.test()
#  PIE.graphics <- function(eval.df, response, grouping.factors, seed=314) {
#    require(dplyr); require(tidyr); require(ggplot2); require(grid); require(gridExtra);
#    ####Data prep
#    set.seed(seed)
#    eval$random <- runif(nrow(eval))
#    Eval <- eval %>%
#      mutate(
#        random = row_number(random)
#      ) %>%
#      gather(var, value, starts_with("boot"), starts_with("KF")) %>%
#      mutate(
#        simMethod = sub("[[:punct:]][0-9A-Za-z]*", "", var),
#        stat = sub("[0-9A-Za-z]*[[:punct:]]", "", var)
#      ) %>%
#      select(-var) %>%
#      spread(stat, value) %>%
#      group_by(.newID) %>%
#      arrange(simMethod) %>%
#      mutate(
#        x=row_number(simMethod),
#        x=(x-mean(x))/10
#      ) %>%
#      group_by(simMethod) %>%
#      mutate(
#        x=row_number(random)+x
#      )
#
#    ####Direct Comparison Plot
#    Eval.small <- arrange(Eval, random) %>% filter(random<=30) %>% arrange(x)
#    p1 <- ggplot(aes(y=fit,x=x, color=simMethod), data=Eval.small) +
#            geom_point(size=I(3)) +
#            geom_linerange(aes(ymax=upr, ymin=lwr), size=I(1)) +
#            geom_point(shape="with.u", color="black", size=I(4),
#                       data=summarize(group_by(Eval.small, .newID), x=mean(x), fit=mean(with.u))) +
#            geom_point(shape="no.u", color="black", size=I(4),
#                       data=summarize(group_by(Eval.small, .newID), x=mean(x), fit=mean(no.u))) +
#            theme_bw() +
#            labs(x="Index", y="Prediction Interval", title="95% Prediction interval by method for 30 random obs") +
#            scale_x_discrete(breaks=1:30, labels=Eval.small$.newID[Eval.small$simMethod=="KF"]) +
#            theme(axis.text.x = element_text(angle=90))
#
#    #####Distribution of fitted values
#    mean.statment <- paste("mean(",response,")", sep="")
#    p2 <- ggplot(aes(x=fit), data=Eval) +
#            geom_density(aes(color=simMethod)) +
#            geom_vline(x=mean(unlist(Eval[,response]))) +
#            geom_density(data=summarize_(
#                                group_by_(Eval, grouping.factors),
#                                fit = mean.statment), color="black") +
#            theme_bw() +
#            labs(x = "Estimate", y="Density",
#                 title="Average point estimates by prediction type \n (vertical black line is sample grand mean, \n black density is distribution of sample group means)")
#
#    #####Bar Graph of KF PI containing other point estimates
#    p3 <- eval %>% select(starts_with("contain.")) %>%
#      gather(simMethod, value) %>%
#      mutate(
#        simMethod = sub("[[:alpha:]]*.", "", simMethod)
#      ) %>%
#      group_by(simMethod) %>%
#      summarize(
#        value=100*mean(value)
#      ) %>%
#      qplot(value, x=simMethod, fill=simMethod, data=., geom="bar", position="dodge", stat="identity") +
#        labs(y="Percent", title="Percent of observations in which our P.I. contains other point estimates") +
#        theme_bw()
#
#    ##Summarizing Bias
#    SD <- sd(eval[,response])
#    bias.data <- eval %>% select(starts_with("distance.")) %>%
#      gather(simMethod, value) %>%
#      mutate(
#        simMethod = sub("[[:alpha:]]*.", "", simMethod),
#        value = value/SD,
#        ymax = max(density(value)$y),
#        xmin = min(value)
#      ) %>%
#      group_by(simMethod) %>%
#      summarize(
#        mean.Distance = round(mean(value),4),
#        mad.Distance = round(mad(value),2),
#        ymax = max(ymax),
#        xmin = min(xmin)
#      )
#    bias.xmin <- min(bias.data$xmin, na.rm=TRUE)
#    bias.xmax <- 0.2 * bias.xmin
#    bias.ymax <- max(bias.data$ymax, na.rm=TRUE)
#    bias.ymin <- .5 * bias.ymax
#    bias.data <- bias.data[,1:3]
#
#    p4a <- eval %>% select(starts_with("distance.")) %>%
#      gather(simMethod, value) %>%
#      mutate(
#        simMethod = sub("[[:alpha:]]*.", "", simMethod),
#        value = value/SD
#      ) %>%
#      qplot(x=value, color=simMethod, data=., geom="density") +
#        labs(x="Standard deviations of response variable") +
#        theme_bw()
#
#    p4b <- tableGrob(bias.data, show.rownames = FALSE)
#    p4 <- arrangeGrob(p4a, p4b, ncol=1, main="Distribution of distance from KF point estimates to other point estimates")
#
#
#    ##Summarizing Bias
#    coverage.data <- eval %>%
#      select(ends_with(".coverage")) %>%
#      gather(simMethod, value) %>%
#      mutate(
#        simMethod = gsub(".coverage","", simMethod, fixed=TRUE),
#        TotalCoverage=value==1,
#        Cov.90to100 = value>0.9 & value < 1,
#        Cov.80to90 = value>0.8 & value <= 0.9,
#        Cov.50to80 = value>0.5 & value <= 0.8,
#        Cov.0to50 = value>0 & value <= 0.5,
#        ZeroCoverage= value==0
#      ) %>%
#      group_by(simMethod) %>%
#      summarize(
#        TotalCoverage= round(100*sum(TotalCoverage)/n(),1),
#        Cover.90to100 =  round(100*sum(Cov.90to100)/n(),1),
#        Cover.80to90 = round(100*sum(Cov.80to90)/n(),1),
#        Cover.50to80 = round(100*sum(Cov.50to80)/n(),1),
#        Cover.0to50 = round(100*sum(Cov.0to50)/n(),1),
#        ZeroCoverage = round(100*sum(ZeroCoverage)/n(),1)
#      )
#    p5a <- eval %>%
#      select(ends_with(".coverage")) %>%
#      gather(simMethod, value) %>%
#      mutate(
#        simMethod = gsub(".coverage","", simMethod, fixed=TRUE)
#      ) %>%
#      ggplot(data=., aes(x=100*value, fill=simMethod)) +
#        geom_bar(aes(y=300*(..count..)/sum(..count..)), binwidth=5) +
#        facet_wrap(~simMethod, ncol=3) +
#        labs(y="Percent of Observations", x="Coverage Percentage") +
#        theme_bw()
#    p5b <- tableGrob(coverage.data, show.rownames = FALSE)
#    p5 <- arrangeGrob(p5a, p5b, ncol=1, main="Distribution of PI coverage percentages")
#
#    ##Wrap-up
#    return(
#      list(
#        CompareRandomObs=p1,
#        FitDistributions=p2,
#        PointEstimateCoverage=p3,
#        BiasSummary=p4,
#        CoverageSummary=p5
#        )
#      )
#  }
# # PIE.graphics(cannonical.1$evalData, response="Reaction", grouping.factors = "Subject")
#
#
# #Step 3: Summarize and Compare Results####
#  #debug(predictInterval.test)
#  cannonical.1 <-  predictInterval.test(model.form = m1.form, model.df = m1.df, predict.df = m1.new.df, nSims=2500)
#  pbPost("note", title="Finished Cannonical Model 1", body=paste(kable(cannonical.1$compareTimes[,1:3]), collapse ="\n"))
#  cannonical.2 <-  predictInterval.test(model.form = m2.form, model.df = m2.df, predict.df = m2.new.df, nSims=2500)
#  pbPost("note", title="Finished Cannonical Model 2", body=paste(kable(cannonical.2$compareTimes[,1:3]), collapse ="\n"))
#  cannonical.4 <-  predictInterval.test(model.form = m4.form, model.df = m4.df, predict.df = m4.new.df, nSims=2500)
#  pbPost("note", title="Finished Cannonical Model 4", body=paste(kable(cannonical.4$compareTimes[,1:3]), collapse ="\n"))
#  cannonical.3 <-  predictInterval.test(model.form = m3.form, model.df = m3.df, model.type="glmer",
#                                        predict.df=m3.new.df, predict.type="response", nSims=2500)
#  pbPost("note", title="Finished Cannonical Model 3", body=paste(kable(cannonical.3$compareTimes[,1:3]), collapse ="\n"))
#
#
#  # save.image()
#
# #Step 4: checking math with a fine-toothed comb ####
#  ##Pull pieces to create toy example
#  KF <- cannonical.1$kf.method
#  model <- cannonical.1$model
#  nsim=5
