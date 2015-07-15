#' Predict from merMod objects with a prediction interval
#'
#' @param model a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param level the width of the prediction interval
#' @param n.sims number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param type type of prediction to develop
#' @param include.resid.var logical, include or exclude the residual varaince for
#' linear models
#' @return 'newdata' with three columns appended, stat and the lower
#'         and upper prediction interval boundaries
#' @export
#' @details Sampling strategy description.
#' @importFrom mvtnorm rmvnorm
#' @import lme4
#' @importFrom abind abind
predictInterval <- function(model, newdata, level = 0.95,
                            n.sims=100, stat=c("median","mean"),
                            type=c("linear.prediction", "probability"),
                            include.resid.var=TRUE){
  outs <- newdata
  predict.type <- match.arg(type,
                            c("linear.prediction", "probability"),
                            several.ok = FALSE)
  stat.type <- match.arg(stat,
                         c("median","mean"),
                         several.ok = FALSE)

  ##First: check if it is a GLMM or NLMM and draw from sigma distribution or incorporate scale parameter if GLMM
  model.devcomp <- getME(model, "devcomp")
  if (model.devcomp$dims[["GLMM"]] == 0 &
      model.devcomp$dims[["NLMM"]] == 0) {
    sigmahat <-  sqrt(1/rgamma(n.sims, 0.5*residDF.merMod(model), 0.5*model.devcomp$cmp[["pwrss"]]))
    if (predict.type=="probability") {
      predict.type="linear.prediction"
      warning("    Asking for predictions on the probability scale makes no sense, resetting predict.type to linear.prediction",
              call.=FALSE)
    }
  }
  else if (model.devcomp$dims[["GLMM"]] == TRUE &
           model@resp$family$family == "binomial" &
           model@resp$family$link == "logit") {
    sigmahat <- rep(1,n.sims)
  }
  else {
    stop("    Prediction for this NLMMs or GLMMs that are not mixed logistic regressions is not yet implemented.")
  }

  # If we do it this way, the function will fail on new data without all levels
  # of X
  # newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X
  # To be sensitive to this, we can take a performance hit and do:
  if(identical(newdata, model@frame)){
    newdata.modelMatrix <- model.matrix(nobars(model@call$formula),
                                        data = model@frame)
  } else{
    tmp <- plyr::rbind.fill(newdata, trimModelFrame(model@frame))
    # attempt to make insensitive to spurious factor levels in betas
    #     nums <- sapply(data, is.numeric); vars <- names(nums[!nums == TRUE])
    #     tmp[, vars] <- apply(tmp[, vars], 2, as.character)
    newdata.modelMatrix <- model.matrix(nobars(model@call$formula),
                                        data = tmp)[1:nrow(newdata), ]
    rm(tmp)
  }

  ##This chunk of code draws from empirical bayes distributions for each level of the random effects
  ##and merges it back onto newdata.
  ##
  ##Right now I am not multiplying the BLUP variance covariance matrices by our
  ##draw of sigma (for linear models) because their variation is unique.  If anything,
  ##this is where one would multiply them by draws of theta from the model.
  reTerms <- reTermNames(model)
  n.reTerms = reTermCount(model)
  reSimA <- array(data = NA, dim = c(nrow(newdata), n.sims, n.reTerms))
  reCoef <- vector(getME(model, "n_rfacs"), mode = "list")
  names(reCoef) <- names(ngrps(model))

  for(j in names(reCoef)){
    reMeans <- as.matrix(ranef(model)[[j]])
    reMatrix <- attr(ranef(model, condVar=TRUE)[[j]], which = "postVar")
    reSimA <- array(data = NA, dim = c(nrow(reMeans), ncol(reMeans), n.sims))
    for(k in 1:nrow(reMeans)){
      meanTmp <- as.matrix(reMeans[k, ])
      matrixTmp <- as.matrix(reMatrix[,,k])
      reSimA[k, ,] <- mvtnorm::rmvnorm(n.sims,
                                       mean=meanTmp,
                                       sigma=matrixTmp,
                                       method="svd")
      dimnames(reSimA)[[1]] <- rownames(ranef(model)[[j]])
      dimnames(reSimA)[[2]] <- colnames(ranef(model)[[j]])
    }
#     reCoef[[j]] <- apply(reSimA, 3,
#                          function(x) merge(newdata, as.data.frame(x),
#                                            by.x = j, by.y = "row.names"))
     reCoef[[j]] <- reSimA
  }

# multiple each simulation of the array by the newdata.ModelMatrix using rbind
  # for missing columns
  # reCoef is now a list that contains all the necessary simulations
  re.xb <- vector(length = length(names(reCoef)), mode = "list")
  names(re.xb) <- names(reCoef)
  for(j in names(reCoef)){
    tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
    keep <- names(tmp)[names(tmp) %in% dimnames(reCoef[[j]])[[2]]]
    tmp <- tmp[, c(keep, "var")]
    tmp$var <- as.character(tmp$var)
    tmpCoef <- reCoef[[j]][, keep, , drop = FALSE]

    tmp.pred <- function(data, coefs){
      yhatTmp <- array(data = NA, dim = c(nrow(data), ncol(data)-1, dim(coefs)[3]))
      for(k in 1:ncol(data)-1){
        for(i in 1:nrow(data)){
          lvl <- as.character(data[, "var"][i])
          yhatTmp[i, k,] <- as.numeric(data[i, k]) * as.numeric(coefs[lvl, k, ])
        }
      }
      return(yhatTmp)
    }
    re.xb[[j]] <- apply(tmp.pred(data = tmp, coefs = tmpCoef), c(1,3), sum)
  }

  ##This chunk of code produces matrix of linear predictors created from the fixed coefs
  ##and incorporate the model's residual variation if requested
  if (include.resid.var==FALSE) {
    if (length(new.levels)==0)
      sigmahat <- rep(1,n.sims)
    else {
      include.resid.var=TRUE
      warning("    \n  Since new levels were detected resetting include.resid.var to TRUE.")
    }
  }

  # fixed.xb is nrow(newdata) x n.sims
  ##Calculate yhat as sum of the components (fixed plus all groupling factors)
  # apply(reSim, c(1,2), function(x), sum(x,na.rm=TRUE))
  betaSim <- abind(lapply(1:n.sims, function(x) mvtnorm::rmvnorm(1, mean = fixef(model), sigma = sigmahat[x]*as.matrix(vcov(model)))), along=1)
  re.xb$fixed <- newdata.modelMatrix %*% t(betaSim)
  yhat <- Reduce('+', re.xb)

  if (include.resid.var==TRUE)
    yhat <- abind::abind(lapply(1:n.sims, function(x) rnorm(nrow(newdata), yhat[,x], sigmahat[x])), along = 2)

  #Output prediction intervals
  if (stat.type == "median") {
    outs$fit <- apply(yhat,1,function(x) as.numeric(quantile(x, .5, na.rm=TRUE)))
  }
  if (stat.type == "mean") {
    outs$fit <- apply(yhat,1,function(x) mean(x, na.rm=TRUE))
  }
  outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2), na.rm=TRUE)))
  outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2), na.rm=TRUE)))
  if (predict.type == "probability") {
    outs$fit <- model@resp$family$linkinv(outs$fit)
    outs$upr <- model@resp$family$linkinv(outs$upr)
    outs$lwr <- model@resp$family$linkinv(outs$lwr)
  }
  #Close it out
  return(outs[, c("fit", "lwr", "upr")])
}
