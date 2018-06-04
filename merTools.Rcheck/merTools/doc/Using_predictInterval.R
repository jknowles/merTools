## ----setup, echo = FALSE, message=FALSE, warning=FALSE, results='hide', cache=FALSE----
library(ggplot2); library(knitr); library(merTools)
knitr::opts_chunk$set(
  cache=FALSE,
  comment="#>",
  collapse=TRUE,
  echo = TRUE
)

## ----Prep, message=FALSE, warning=FALSE----------------------------------
set.seed(271828)
data(sleepstudy)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
display(fm1)

## ----predInt-------------------------------------------------------------
PI.time <- system.time(
  PI <- predictInterval(merMod = fm1, newdata = sleepstudy,
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
)

## ----Inspect predInt, results="asis", echo=FALSE-------------------------
kable(head(PI))

## ----Inspect predInt 2, fig.width=7, fig.align="center"------------------
library(ggplot2);
ggplot(aes(x=1:30, y=fit, ymin=lwr, ymax=upr), data=PI[1:30,]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") + theme_bw()

## ----PrepFat, message=FALSE, warning=FALSE-------------------------------
ais <- DAAG::ais
fatmodel = lmer(pcBfat~ ht + sex + (1 | sport),data=ais)

## ----ShowFat, message=FALSE, warning=FALSE-------------------------------
ais[1,]
predictInterval(fatmodel, ais[1,], include.resid.var=0) #predict the average body fat for a group of 196cm female baseball players

## ----ShowFat2, message=FALSE, warning=FALSE------------------------------
predictInterval(fatmodel, ais[1,], include.resid.var=0, ignore.fixed.terms = 1) #predict the average body fat for a group of 196cm female baseball players, taking the global intercept (mean body fat) as fully known
predictInterval(fatmodel, ais[1,], include.resid.var=0, ignore.fixed.terms = "(Intercept)") #Same as above
predictInterval(fatmodel, ais[1,], include.resid.var=0, ignore.fixed.terms = 1:2) #as above, taking the first two fixed effects (intercept and height effect) as fully known

## ----ShowFat3, message=FALSE, warning=FALSE------------------------------
predictInterval(fatmodel, ais[1,], include.resid.var=0, fix.intercept.variance = T) #predict the average body fat for a group of 196cm female baseball players, using an ad-hoc correction for the covariance of the intercept with the random intercept effects.

## ----arm.Sim, fig.width=7, fig.height=4, fig.align="center"--------------
PI.arm.time <- system.time(
  PI.arm.sims <- arm::sim(fm1, 1000)
)

PI.arm <- data.frame(
  fit=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.500)),
  upr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.975)),
  lwr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.025))
)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="arm::sim()", x=(1:nrow(PI.arm))+0.1, PI.arm))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)



## ----bootMer.1, fig.width=7, fig.height=4, fig.align="center"------------
##Functions for bootMer() and objects
####Return predicted values from bootstrap
mySumm <- function(.) {
  predict(., newdata=sleepstudy, re.form=NULL)
}
####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

##lme4::bootMer() method 1
PI.boot1.time <- system.time(
  boot1 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=FALSE, type="parametric")
)

PI.boot1 <- sumBoot(boot1)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 1", x=(1:nrow(PI.boot1))+0.1, PI.boot1))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)


## ----bootMer.2, fig.width=7, fig.height=4, fig.align="center"------------
##lme4::bootMer() method 2
PI.boot2.time <- system.time(
  boot2 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=TRUE, type="parametric")
)

PI.boot2 <- sumBoot(boot2)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 2", x=(1:nrow(PI.boot2))+0.1, PI.boot2))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)



## ----bootMer.3, fig.width=7, fig.height=4, fig.align="center"------------
##lme4::bootMer() method 3
PI.boot3.time <- system.time(
  boot3 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=TRUE, type="semiparametric")
)

PI.boot3 <- sumBoot(boot3)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 3", x=(1:nrow(PI.boot3))+0.1, PI.boot3))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)



## ----echo=FALSE, message=FALSE-------------------------------------------
library(rstanarm)

central_intervals <- function(x, prob) {
  if (!identical(length(prob), 1L) || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.",
         call. = FALSE)
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  out <- t(apply(x, 2L, quantile, probs = probs))
  structure(out, dimnames = list(colnames(x), labs))
}


## ---- message=FALSE, fig.width=7, fig.height=4, fig.align="center"-------
PI.time.stan <- system.time({
  fm_stan <- stan_lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy,
                       verbose = FALSE, open_progress = FALSE, refresh = -1,
                       show_messages=FALSE)
  zed <- posterior_predict(fm_stan)
  PI.stan <- cbind(apply(zed, 2, median), central_intervals(zed, prob=0.95))
})


print(fm_stan)

PI.stan <- as.data.frame(PI.stan)
names(PI.stan) <- c("fit", "lwr", "upr")
PI.stan <- PI.stan[, c("fit", "upr", "lwr")]
comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="rstanArm", x=(1:nrow(PI.stan))+0.1, PI.stan))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() +
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)


## ---- echo=FALSE---------------------------------------------------------
times <- rbind(PI.time, PI.arm.time, PI.boot1.time, PI.boot2.time, PI.boot3.time,
               PI.time.stan)[, 1:3]
rownames(times) <- c("predictInterval()", "arm::sim()", "lme4::bootMer()-Method 1", "lme4::bootMer()-Method 2", "lme4::bootMer()-Method 3", "rstanarm:predict")
kable(times)

