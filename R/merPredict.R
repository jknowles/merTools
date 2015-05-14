#' Predict from merMod objects with a prediction interval
#'
#' @param model a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param level the width of the prediction interval
#' @param nsim number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param predict.type type of prediction to develop
#'
#' @return The original data with three columns of predictors appended
#' @export
#' @details Sampling strategy description.
#' @importFrom mvtnorm rmvnorm
#' @import lme4
predictInterval <- function(model, newdata, level = 0.95,
                            nsim=1000, stat="median",
                            predict.type=c("link", "response")){
  outs <- newdata
  predict.type <- match.arg("predict.type",
                            c("link", "response"),
                            several.ok = FALSE)
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
    reSim[[group]] <- merge(newdata, reSim[[group]], by=group, all.x=TRUE)
    reSim[[group]] <- as.matrix(reSim[[group]][,setdiff(cnames, group)])
  }

  #Calculate yhat as sum of components
  yhat <- fixed.xb + apply(simplify2array(reSim), c(1,2), sum)

  #Output prediction intervals
  if (stat == "median") {
    outs$fit <- apply(yhat,1,function(x) as.numeric(quantile(x, .5)))
  }
  if (stat == "mean") {
    outs$fit <- apply(yhat,1,mean)
  }
  outs$upr <- apply(yhat,1,function(x) as.numeric(quantile(x, 1 - ((1-level)/2))))
  outs$lwr <- apply(yhat,1,function(x) as.numeric(quantile(x, ((1-level)/2))))
  if (predict.type == "response") {
    outs$fit <- model@resp$family$linkinv(outs$fit)
    outs$upr <- model@resp$family$linkinv(outs$upr)
    outr$lwr <- model@resp$family$linkinc(outs$lwr)
  }
  #Close it out
  return(outs)
}
