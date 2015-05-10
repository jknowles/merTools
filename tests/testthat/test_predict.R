# # -----------------------------------------------------
# #-------------------------------------------------------
# set.seed(51315)
# library(lme4)
# data(grouseticks)
# grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
# grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
#
#
# # Build out models
# form <- TICKS ~ YEAR + HEIGHT +(1|BROOD) + (1|LOCATION) + (1|INDEX)
# glmer3Lev  <- glmer(form, family="poisson",data=grouseticks,
#                     control = glmerControl(optimizer="Nelder_Mead",
#                                            optCtrl=list(maxfun = 1e5)))
#
# # Sleepstudy
# lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#
# # Wackier example
# data(Orthodont,package="nlme")
# Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
# Orthodont$nsexage <- with(Orthodont, nsex*age)
# lmerSlope2 <- lmer(distance ~ age + (0 + age + nsex|Subject), data=Orthodont)
#
# uj <- REsim(lmerSlope1, nsims = 1000)
# betas <- FEsim(lmerSlope1, nsims = 1000)
# newdata <- sleepstudy[1:100,]
#
# library(mvtnorm)

#
# predictInterval <- function(model, newdata, level = 0.95){
#   require(mvtnorm)
#   reTerms <- names(ngrps(model))
#   if(length(reTerms) > 1){
#     stop("Multiple grouping terms not yet implemented.")
#   }
#   group <- reTerms[1]
#   outs <- data.frame(fit = rep(NA, nrow(newdata)),
#                      lwr = rep(NA, nrow(newdata)),
#                      upr = rep(NA, nrow(newdata)))
#
#   reMatrixExtract <- function(model, group, case){
#     slotNum <- ranef(model)[[group]]
#     slotNum$seq <- 1:nrow(slotNum)
#     slotNum <- slotNum[case, ]$seq
#     mat <- as.matrix(attr(ranef(model, condVar=TRUE)[[group]],
#                           which = "postVar")[, , slotNum])
#     row.names(mat) <- names(ranef(model)[[group]])
#     colnames(mat) <- row.names(mat)
#     return(mat)
#   }
#
#   reMeansExtract <- function(model, group, case){
#     rerow <- ranef(model)[[group]][case, ]
#     out <- as.numeric(rerow)
#     names(out) <- names(rerow)
#     return(out)
#   }
#   combineREFE <- function(reSim, betaSim){
#     comboCols <- union(colnames(reSim), colnames(betaSim))
#     REonly <- setdiff(colnames(reSim), colnames(betaSim))
#     FEonly <- setdiff(colnames(betaSim), colnames(reSim))
#     totCoef <- reSim[, comboCols] + betaSim[, comboCols]
#     totCoef <- cbind(totCoef, reSim[, REonly])
#     totCoef <- cbind(totCoef, betaSim[, FEonly])
#     return(totCoef)
#   }
#
#   for(i in 1:nrow(newdata)){
#     caseID <- newdata[i, group]
#     betaSim <- rmvnorm(1000, mean = fixef(model),
#                        sigma = as.matrix(vcov(model)))
#     reSim <- rmvnorm(1000, mean = reMeansExtract(model, group = group, case=caseID),
#                      sigma = reMatrixExtract(model, group = group, case = caseID))
#     totCoef <- combineREFE(reSim, betaSim)
#     yhat <- newdata[i, "Days"] * totCoef[, "Days"] + totCoef[, "(Intercept)"]
#     outs[i, "fit"] <- mean(yhat)
#     outs[i, "upr"] <- as.numeric(quantile(yhat, level))
#     outs[i, "lwr"] <- as.numeric(quantile(yhat, 1 - level))
#   }
#   return(outs)
# }
#
# predExamp <- cbind(predictInterval(lmerSlope1, sleepstudy[1:10, ], level = 0.95),
#                    predict(lmerSlope1, sleepstudy[1:10,]),
#                    predict(lmerSlope1, sleepstudy[1:10,], re.form = NA))







