# substEffSimFE <- function(mod, var, ncases, unscale = NULL, unscaleDV = NULL, ...){
#   if(missing(unscaleDV)){
#     if(class(mod) == "lmerMod"){
#       unscaleDV <- TRUE
#     } else {
#       unscaleDV <- FALSE
#     }
#   }
#   dvName <- names(mod@frame)[1]
#   # var should be a character, name of variable in model
#   # ncases is integer
#   # unscale is a logical indicating if the independent variable should be unscaled
#   # mod is a merMod
#   cases <- mod@frame[sample(1:nrow(mod@frame), ncases),]
#   cases$caseID <- 1:nrow(cases)
#   varLength <- length(unique(mod@frame[, var]))
#   if(varLength < 30){
#     jitters <- unique(mod@frame[, var])
#   } else {
#     jitters <- quantile(mod@frame[, var], probs=seq(0,1, by =0.05))
#   }
#   simvals <- expand.grid(caseID = unique(cases$caseID), newvar = jitters)
#   plotdf <- merge(cases, simvals)
#   plotdf$oldvar <- plotdf[, var]
#   plotdf[, var] <- plotdf$newvar
#   plotdf$newvar <- NULL
#   plotdf <- cbind(plotdf, easyPredCI(mod, newdata=plotdf, ...))
#   if(missing(unscale)){
#     return(plotdf)
#   } else {
#     stopifnot(class(unscale) == "matrix")
#     plotdf$newvarUnscale <- (plotdf[, var] * 2 * unscale[var, 2]) + unscale[var, 1]
#     if(unscaleDV == TRUE){
#       plotdf$yhatUnscale <- (plotdf[, "yhat"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
#       plotdf$lwrUnscale <- (plotdf[, "lwr"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
#       plotdf$uprUnscale <- (plotdf[, "upr"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
#     }
#     return(plotdf)
#   }
# }
#
#
# substEffSimRE <- function(mod, var, ncases, ...){
#   # var should be a character, name of variable in model
#   # ncases is integer
#   # mod is a merMod
#   cases <- mod@frame[sample(1:nrow(mod@frame), ncases),]
#   cases$caseID <- 1:nrow(cases)
#   varLength <- length(unique(mod@frame[, var]))
#   jitters <- unique(mod@frame[, var])
#   simvals <- expand.grid(caseID = unique(cases$caseID), newvar = jitters)
#   plotdf <- merge(cases, simvals)
#   plotdf$oldvar <- plotdf[, var]
#   plotdf[, var] <- plotdf$newvar
#   plotdf <- cbind(plotdf, easyPredCI(mod, newdata=plotdf, re = TRUE, ...))
#   return(plotdf)
# }
#
#
# easyPredCI <- function(model, newdata, alpha=0.05, re = NULL) {
#   # From https://rpubs.com/bbolker/glmmchapter
#   ## baseline prediction, on the linear predictor (logit) scale:
#   if(missing(re)){
#     pred0 <- predict(model, re.form = NA, newdata=newdata)
#   } else{
#     pred0 <- predict(model, newdata=newdata)
#   }
#   ## fixed-effects model matrix for new data
#   X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
#                     newdata)
#   beta <- fixef(model) ## fixed-effects coefficients
#   V <- vcov(model)     ## variance-covariance matrix of beta
#   pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
#
#   if(class(model) == "lmerMod"){
#     linkinv <- identity
#   } else if(class(model) == "glmerMod"){
#     ## inverse-link (logistic) function: could also use plogis()
#     linkinv <- model@resp$family$linkinv
#   }
#   ## construct 95% Normal CIs on the link scale and
#   ##  transform back to the response (probability) scale:
#   crit <- -qnorm(alpha/2)
#   linkinv(cbind(lwr = pred0 - crit * pred.se,
#                 upr = pred0 + crit * pred.se,
#                 yhat = pred0))
# }
#
