#' @title Extracts random effects
#' @name REextract
#' @description Extracts random effect terms from an lme4 model
#' @param mod a merMod object from the lme4 package
#' @importFrom plyr adply
#' @return a data frame with random effects from lme4
#' @export
REextract <- function(mod){
  out <- lme4::ranef(mod, condVar = TRUE)
  out.se <- plyr::adply(attr(out[[1]], which = "postVar"), c(3),
                        function(x) sqrt(diag(x)))
  out.pt <- out[[1]]
  names(out.se)[-1] <- paste0(names(out.pt), "_se")
  newout <- cbind(out.pt, out.se[, -1])
  return(newout)
}

#' @title Simulate random effects from merMod
#' @name REsim
#' @description Simulate random effects from merMod object posterior distributions
#' @param x a merMod object from the lme4 package
#' @param nsims number of simulations to use
#' @param OR logical, should parameters be converted to odds ratios?
#' @importFrom arm sim
#' @import lme4
#' @return a data frame with distribution of random effect parameters
#' @details Use the Gelman sim technique to build empirical Bayes estimates.
#'  Uses the sim function in the arm package
#' @export
REsim <- function(x, nsims, OR = FALSE){
  mysim <- sim(x, n.sims = nsims)
  reDims <- length(mysim@ranef)
  tmp.out <- vector("list", reDims)
  for(i in c(1:reDims)){
    tmp.out[[i]] <- plyr::adply(mysim@ranef[[i]], c(2, 3), plyr::each(c(mean, median, sd)))
    tmp.out[[i]]$level <- paste0("Level ", i)
    tmp.out[[i]]$level <- as.character(tmp.out[[i]]$level)
    tmp.out[[i]]$X1 <- as.character(tmp.out[[i]]$X1)
    tmp.out[[i]]$X2 <- as.character(tmp.out[[i]]$X2)
  }
  dat <- do.call(rbind, tmp.out)
  if(OR == TRUE){
    dat$median <- exp(dat$median)
    dat$mean <- exp(dat$mean)
    dat$sd <- NA # don't know how to do SE of odds ratios currently
    return(dat)
  } else{
    return(dat)
  }
}
#
#
# FEsim <- function(x, nsims){
#   mysim <- sim(x, n.sims = nsims)
#   means <- apply(mysim@fixef, MARGIN = 2, mean)
#   medians <- apply(mysim@fixef, MARGIN = 2, median)
#   sds <- apply(mysim@fixef, MARGIN =2, sd)
#   dat <- data.frame(var = names(means), meanEff = means, medEff = medians,
#                     sdEff = sds, row.names=NULL)
#   return(dat)
# }
#
# # http://stats.stackexchange.com/questions/56525/standard-deviation-of-random-effect-is-0
# # http://stat.columbia.edu/~jcliu/paper/HierarchicalPrior.pdf
# # http://www.stat.columbia.edu/~radon/paper/paper.pdf (example)
#
# # Plot the effects
# plotMCMCre <- function(dat, scale, var, sd, sigmaScale = NULL, oddsRatio = FALSE,
#                        labs = NULL){
#   if(!missing(sigmaScale)){
#     dat[, sd] <- dat[, sd] / sigmaScale
#     dat[, var] <- dat[, var] / sigmaScale
#   }
#   dat[, sd] <- dat[, sd] * scale
#   dat[, "ymax"] <- dat[, var] + dat[, sd]
#   dat[, "ymin"] <- dat[, var] - dat[, sd]
#   hlineInt <- 0
#   if(oddsRatio == TRUE){
#     dat[, "ymax"] <- exp(dat[, "ymax"])
#     dat[, var] <- exp(dat[, var])
#     dat[, "ymin"] <- exp(dat[, "ymin"])
#     hlineInt <- 1
#   }
#   if(!missing(labs)){
#     xlabs.tmp <- element_text(face = "bold")
#     xvar <- labs
#   } else {
#     xlabs.tmp <- element_blank()
#     xvar <- "id"
#   }
#
#   dat[order(dat[, var]), "id"] <- c(1:nrow(dat))
#   ggplot(dat, aes_string(x = xvar, y = var, ymax = "ymax",
#                          ymin = "ymin")) +
#     geom_pointrange(alpha = I(0.4)) + theme_dpi() + geom_point() +
#     labs(x = "Group", y = "Effect Range", title = "Effect Ranges") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.text.x = xlabs.tmp, axis.ticks.x = element_blank()) +
#     geom_hline(yintercept = hlineInt, color = I("red"), size = I(1.1))
# }

plotMCMCre <- function(dat, scale, var, sd, sigmaScale = NULL, oddsRatio = FALSE, 
                       labs = NULL){
  if(!missing(sigmaScale)){
    dat[, sd] <- dat[, sd] / sigmaScale
    dat[, var] <- dat[, var] / sigmaScale
  }
  dat[, sd] <- dat[, sd] * scale
  dat[, "ymax"] <- dat[, var] + dat[, sd] 
  dat[, "ymin"] <- dat[, var] - dat[, sd]
  hlineInt <- 0
  if(oddsRatio == TRUE){
    dat[, "ymax"] <- exp(dat[, "ymax"])
    dat[, var] <- exp(dat[, var])
    dat[, "ymin"] <- exp(dat[, "ymin"])
    hlineInt <- 1
  }
  if(!missing(labs)){
    xlabs.tmp <- element_text(face = "bold")
    xvar <- labs
  } else {
    xlabs.tmp <- element_blank()
    xvar <- "id"
  }
  
  dat[order(dat[, var]), "id"] <- c(1:nrow(dat))
  ggplot(dat, aes_string(x = xvar, y = var, ymax = "ymax", 
                         ymin = "ymin")) + 
    geom_pointrange(alpha = I(0.4)) + theme_dpi() + geom_point() +
    labs(x = "Group", y = "Effect Range", title = "Effect Ranges") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = xlabs.tmp, axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = hlineInt, color = I("red"), size = I(1.1))
}


RMSE.loess <- function(m){
  sqrt(sum(m$residuals^2)/length(m$residuals))
}

substEffSimFE <- function(mod, var, ncases, unscale = NULL, unscaleDV = NULL, ...){
  if(missing(unscaleDV)){
    if(class(mod) == "lmerMod"){
      unscaleDV <- TRUE
    } else {
      unscaleDV <- FALSE
    }
  }
  dvName <- names(mod@frame)[1]
  # var should be a character, name of variable in model
  # ncases is integer
  # unscale is a logical indicating if the independent variable should be unscaled
  # mod is a merMod
  cases <- mod@frame[sample(1:nrow(mod@frame), ncases),]
  cases$caseID <- 1:nrow(cases)
  varLength <- length(unique(mod@frame[, var]))
  if(varLength < 30){
    jitters <- unique(mod@frame[, var])
  } else {
    jitters <- quantile(mod@frame[, var], probs=seq(0,1, by =0.05))
  }
  simvals <- expand.grid(caseID = unique(cases$caseID), newvar = jitters)
  plotdf <- merge(cases, simvals)
  plotdf$oldvar <- plotdf[, var]
  plotdf[, var] <- plotdf$newvar
  plotdf$newvar <- NULL
  plotdf <- cbind(plotdf, easyPredCI(mod, newdata=plotdf, ...))
  if(missing(unscale)){
    return(plotdf)
  } else {
    stopifnot(class(unscale) == "matrix")
    plotdf$newvarUnscale <- (plotdf[, var] * 2 * unscale[var, 2]) + unscale[var, 1]
    if(unscaleDV == TRUE){
      plotdf$yhatUnscale <- (plotdf[, "yhat"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
      plotdf$lwrUnscale <- (plotdf[, "lwr"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
      plotdf$uprUnscale <- (plotdf[, "upr"] * 2 * unscale[dvName, 2]) + unscale[dvName, 1]
    }
    return(plotdf)
  }
}


substEffSimRE <- function(mod, var, ncases, ...){
  # var should be a character, name of variable in model
  # ncases is integer
  # mod is a merMod
  cases <- mod@frame[sample(1:nrow(mod@frame), ncases),]
  cases$caseID <- 1:nrow(cases)
  varLength <- length(unique(mod@frame[, var]))
  jitters <- unique(mod@frame[, var])
  simvals <- expand.grid(caseID = unique(cases$caseID), newvar = jitters)
  plotdf <- merge(cases, simvals)
  plotdf$oldvar <- plotdf[, var]
  plotdf[, var] <- plotdf$newvar
  plotdf <- cbind(plotdf, easyPredCI(mod, newdata=plotdf, re = TRUE, ...))
  return(plotdf)
}


easyPredCI <- function(model, newdata, alpha=0.05, re = NULL) {
  # From https://rpubs.com/bbolker/glmmchapter
  ## baseline prediction, on the linear predictor (logit) scale:
  if(missing(re)){
    pred0 <- predict(model, re.form = NA, newdata=newdata)
  } else{
    pred0 <- predict(model, newdata=newdata)
  }
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                    newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  
  if(class(model) == "lmerMod"){
    linkinv <- identity
  } else if(class(model) == "glmerMod"){
    ## inverse-link (logistic) function: could also use plogis()
    linkinv <- model@resp$family$linkinv
  }
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr = pred0 - crit * pred.se,
                upr = pred0 + crit * pred.se, 
                yhat = pred0))
}

RMSE.merMod <- function(mod){
  # Express RMSE as percentage of dependent variable range
  dvSD <- sd(mod@frame[, 1])
  RMSE <- sqrt(mean(residuals(mod)^2))
  return(RMSE/dvSD)
}

bootRMSE <- function(data, indx, model){
  moddat <- data[indx, ]
  modb <- eval(model@call)
  RMSE.merMod(modb)
}
