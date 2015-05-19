#' Predict from merMod objects with a prediction interval
#'
#' @param model a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param level the width of the prediction interval
#' @param nsim number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param predict.type type of prediction to develop
#'
#' @return 'newdata' with three columns appended, stat and the lower
#'         and upper prediction interval boundaries
#' @export
#' @details Sampling strategy description.
#' @importFrom mvtnorm rmvnorm
#' @import lme4
predictInterval <- function(model, newdata, level = 0.95,
                            nsim=1000, stat="median",
                            predict.type=c("link", "response")){
  outs <- newdata
#  predict.type <- match.arg("predict.type",
#                            c("link", "response"),
#                            several.ok = FALSE)
  #Sort out all the levels
  reTerms <- names(ngrps(model))
  n.reTerms = length(reTerms)
  ##The following 3 lines produce a matrix of linear predictors created from the fixed coefs
  betaSim <- rmvnorm(nsim, mean = fixef(model), sigma = as.matrix(vcov(model)))
  newdata.modelMatrix <- lFormula(formula = model@call, data=newdata)$X
  fixed.xb <- newdata.modelMatrix %*% t(betaSim)

  ##Draw from random effects distributions for each level and merge onto newdata
  reSim <- NULL
  for (j in seq_along(reTerms)) {
    group=reTerms[j]
    reMeans <- array(ranef(model)[[group]])
    reMatrix <- attr(ranef(model, condVar=TRUE)[[group]], which = "postVar")
    reSim[[group]] <- data.frame(rownames(reMeans), matrix(NA, nrow=nrow(reMeans),
                                                           ncol=nsim))
    colnames(reSim[[group]]) <- c(group, paste("sim", 1:nsim, sep=""))
    for (k in 1:nrow(reMeans)) {
      lvl = rownames(reMeans)[k]
      reSim[[group]][k,2:ncol(reSim[[group]])] <- rmvnorm(nsim,
                                                          mean=as.matrix(reMeans[k,]),
                                                          sigma=as.matrix(reMatrix[,,k]))
    }
    cnames <- colnames(reSim[[group]])
    #This is where to check for groups that only exist in newdata
    new.levels <- setdiff(newdata[,group], reSim[[group]][,1])
    if (length(new.levels)>0) {
      msg <- paste("     The following levels of ", group, " from newdata -- ", paste0(new.levels, collapse=", "),
                   " -- are not in the model data. \n     Currently, predictions for these values are based only on the fixed coefficients \n     and <<<NOT>>> on the observation-level error.", sep="")
      warning(msg, call.=FALSE)
    }
    reSim[[group]] <- merge(newdata, reSim[[group]], by=group, all.x=TRUE)
    reSim[[group]] <- as.matrix(reSim[[group]][,setdiff(cnames, group)])
  }

  #Calculate yhat as sum of components
  yhat <- fixed.xb + apply(simplify2array(reSim), c(1,2), function(x) sum(x, na.rm=TRUE))

  #Output prediction intervals
  if (stat == "median") {
    outs$fit <- apply(yhat,1,function(x) as.numeric(quantile(x, .5, na.rm=TRUE)))
  }
  if (stat == "mean") {
    outs$fit <- apply(yhat,1,function(x) mean(x, na.rm=TRUE))
  }
  outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2), na.rm=TRUE)))
  outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2), na.rm=TRUE)))
  if (predict.type == "response") {
    outs$fit <- model@resp$family$linkinv(outs$fit)
    outs$upr <- model@resp$family$linkinv(outs$upr)
    outr$lwr <- model@resp$family$linkinc(outs$lwr)
  }
  #Close it out
  return(outs)
}
