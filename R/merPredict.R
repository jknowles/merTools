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
#' @importFrom plyr rbind.fill
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

  # Fix formula to allow for random slopes not in the fixed slopes
  matrixForm <- formulaBuild(model)

  # If we do it this way, the function will fail on new data without all levels
  # of X
  # newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X
  # To be sensitive to this, we can take a performance hit and do:
  if(identical(newdata, model@frame)){
    newdata.modelMatrix <- model.matrix(matrixForm,
                                        data = model@frame)
  } else{
    tmp <- plyr::rbind.fill(newdata, trimModelFrame(model@frame))
    # attempt to make insensitive to spurious factor levels in betas
    #     nums <- sapply(data, is.numeric); vars <- names(nums[!nums == TRUE])
    #     tmp[, vars] <- apply(tmp[, vars], 2, as.character)
    newdata.modelMatrix <- model.matrix(matrixForm,
                                        data = tmp)[1:nrow(newdata), ]
    rm(tmp)
  }

  ##This chunk of code draws from empirical bayes distributions for each level of the random effects
  ##and merges it back onto newdata.
  ##
  ##Right now I am not multiplying the BLUP variance covariance matrices by our
  ##draw of sigma (for linear models) because their variation is unique.  If anything,
  ##this is where one would multiply them by draws of theta from the model.
  re.xb <- vector(getME(model, "n_rfacs"), mode = "list")
  names(re.xb) <- names(ngrps(model))
    for(j in names(re.xb)){
    reMeans <- as.matrix(ranef(model)[[j]])
    reMatrix <- attr(ranef(model, condVar=TRUE)[[j]], which = "postVar")
    # OK, let's knock out all the random effects we don't need
    obslvl <- unique(as.character(newdata[, j]))
    alllvl <- rownames(reMeans)
    keep <- intersect(obslvl, alllvl)
    # Add switch if no random groups are observed to avoid indexing errors,
    # we burn 1 sample of 1 group of all coefficients that will eventually
    # be multiplied by zero later on
    if(length(keep) > 0){
      reMeans <- reMeans[keep, , drop=FALSE]
      dimnames(reMatrix)[[3]] <- alllvl
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else{
      reMeans <- reMeans[1, , drop=FALSE]
      reMatrix <- reMatrix[, , 1, drop = FALSE]
    }
    #
    reSimA <- array(data = NA, dim = c(nrow(reMeans), ncol(reMeans), n.sims))
    for(k in 1:nrow(reMeans)){
      meanTmp <- as.matrix(reMeans[k, ])
      matrixTmp <- as.matrix(reMatrix[,,k])
      reSimA[k, ,] <- mvtnorm::rmvnorm(n.sims,
                                       mean=meanTmp,
                                       sigma=matrixTmp,
                                       method="svd")
      rownames(reSimA) <- rownames(reMeans)
      colnames(reSimA) <- colnames(reMeans)
    }
     tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
     keep <- names(tmp)[names(tmp) %in% dimnames(reSimA)[[2]]]
     tmp <- tmp[, c(keep, "var")]
     tmp$var <- as.character(tmp$var)
     colnames(tmp)[which(names(tmp) == "var")] <- names(newdata[, j,  drop=FALSE])
     tmpCoef <- reSimA[, keep, , drop = FALSE]
     tmp.pred <- function(data, coefs, group){
       yhatTmp <- array(data = NA, dim = c(nrow(data), ncol(data)-1, dim(coefs)[3]))
   new.levels <- unique(as.character(data[, group])[!as.character(data[, group]) %in% dimnames(coefs)[[1]]])
       msg <- paste("     The following levels of ", group, " from newdata \n -- ", paste0(new.levels, collapse=", "),
                    " -- are not in the model data. \n     Currently, predictions for these values are based only on the \n fixed coefficients and the observation-level error.", sep="")
       if(length(new.levels > 0)){
         warning(msg, call.=FALSE)
       }
       for(k in 1:ncol(data)-1){
         for(i in 1:nrow(data)){
           lvl <- as.character(data[, group][i])
           if(lvl %in% dimnames(coefs)[[1]]){
             yhatTmp[i, k,] <- as.numeric(data[i, k]) * as.numeric(coefs[lvl, k, ])
           } else{
             yhatTmp[i, k,] <- as.numeric(data[i, k]) * 0
           }

         }
       }
       return(yhatTmp)
     }
     re.xb[[j]] <- apply(tmp.pred(data = tmp, coefs = tmpCoef,
                                  group = names(newdata[, j, drop=FALSE])), c(1,3), sum)
  }

  # TODO: Add a check for new.levels that is outside of the above loop
  # for now, ignore this check
  if (include.resid.var==FALSE) {
#     if (length(new.levels)==0)
      sigmahat <- rep(1,n.sims)
#     else {
#       include.resid.var=TRUE
#       warning("    \n  Since new levels were detected resetting include.resid.var to TRUE.")
#     }
  }

  # fixed.xb is nrow(newdata) x n.sims
  ##Calculate yhat as sum of the components (fixed plus all groupling factors)
  # apply(reSim, c(1,2), function(x), sum(x,na.rm=TRUE))
  betaSim <- abind::abind(lapply(1:n.sims, function(x) mvtnorm::rmvnorm(1, mean = fixef(model), sigma = sigmahat[x]*as.matrix(vcov(model)))), along=1)
  # Pad betaSim
  if(ncol(newdata.modelMatrix) > ncol(betaSim)){
    pad <- matrix(rep(0), nrow = nrow(betaSim),
                  ncol = ncol(newdata.modelMatrix) - ncol(betaSim))
    if(ncol(pad) > 0){
      message("Fixed effect matrix has been padded with 0 coefficients
            for random slopes not included in the fixed effects and interaction terms.")
    }
    colnames(pad) <- setdiff(colnames(newdata.modelMatrix), colnames(betaSim))
    betaSim <- cbind(betaSim, pad)
    keep <- intersect(colnames(newdata.modelMatrix), colnames(betaSim))
    newdata.modelMatrix <- newdata.modelMatrix[, keep]
    betaSim <- betaSim[, keep]
  }

  re.xb$fixed <- newdata.modelMatrix %*% t(betaSim)
  yhat <- Reduce('+', re.xb)
  # alternative if missing data present:
  # yhat <- apply(simplify2array(re.xb), c(1,2), sum)

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
