#' Predict from merMod objects with a prediction interval
#' @description This function provides a way to capture model uncertainty in
#' predictions from multi-level models fit with \code{lme4}. By drawing a sampling
#' distribution for the random and the fixed effects and then estimating the fitted
#' value across that distribution, it is possible to generate a prediction interval
#' for fitted values that includes all variation in the model except for variation
#' in the covariance paramters, theta. This is a much faster alternative than
#' bootstrapping for models fit to medium to large datasets.
#' @param merMod a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param level the width of the prediction interval
#' @param n.sims number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param type type of prediction to develop
#' @param include.resid.var logical, include or exclude the residual variance for
#' linear models
#' @param returnSims logical, should all n.sims simulations be returned?
#' @param seed numeric, optional argument to set seed for simulations
#' @param .parallel, logical should parallel computation be used, default is FALSE
#' @param .paropts, a list of additional options passed into the foreach function
#' when parallel computation is enabled. This is important if (for example) your
#' code relies on external data or packages: use the .export and .packages arguments
#'  to supply them so that all cluster nodes have the correct environment set up
#'  for computing.
#' @return a data.frame with three columns:
#' \describe{
#'     \item{\code{fit}}{The center of the distribution of predicted values as defined by
#'     the \code{stat} parameter.}
#'     \item{\code{lwr}}{The lower confidence interval bound corresponding to the quantile cut
#'     defined in \code{level}.}
#'     \item{\code{upr}}{The upper confidence interval bound corresponding to the quantile cut
#'     defined in \code{level}.}
#'   }
#' If returnSims = TRUE, then the individual simulations are attached to this
#' data.frame in the attribute \code{sim.results} and are stored as a matrix.
#' @details To generate a prediction inteval, the function first computes a simulated
#' distribution of all of the parameters in the model. For the random, or grouping,
#' effects, this is done by sampling from a multivariate normal distribution which
#' is defined by the BLUP estimate provided by \code{ranef} and the associated
#' variance-covariance matrix for each observed level of each grouping terms. For
#' each grouping term, an array is build that has as many rows as there are levels
#' of the grouping factor, as many columns as there are predictors at that level
#' (e.g. an intercept and slope), and is stacked as high as there are number of
#' simulations. These arrays are then multiplied by the new data provided to the
#' function to produce a matrix of yhat values. The result is a matrix of the simulated
#' values of the linear predictor for each observation for each simulation. Each
#' grouping term has such a matrix for each observation. These values can be added
#' to get the estimate of the fitted value for the random effect terms, and this
#' can then be added to a matrix of simulated values for the fixed effect level to
#' come up with \code{n.sims} number of possible yhat values for each observation.
#'
#' The distribution of simulated values is cut according to the interval requested
#' by the function. The median or mean value as well as the upper and lower bounds
#' are then returned. These can be presented either on the linear predictor scale
#' or on the response scale using the link function in the \code{merMod}.
#' @note \code{merTools} includes the functions \code{subBoot} and \code{thetaExtract}
#' to allow the user to estimate the variability in \code{theta} from a larger
#' model by bootstrapping the model fit on a subset, to allow faster estimation.
#' @export
#' @importFrom mvtnorm rmvnorm
#' @import lme4
#' @import plyr
#' @importFrom abind abind
#' @examples
#' m1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
#' regFit <- predict(m1, newdata = sleepstudy[11, ]) # a single value is returned
#' intFit <- predictInterval(m1, newdata = sleepstudy[11, ]) # bounded values
#' # Can do glmer
#' d1 <- cbpp
#' d1$y <- d1$incidence / d1$size
#'  gm2 <- glmer(y ~ period + (1 | herd), family = binomial, data = d1,
#'                nAGQ = 9, weights = d1$size)
#'  regFit <- predict(gm2, newdata = d1[1:10, ])
#'  # get probabilities
#'  regFit <- predict(gm2, newdata = d1[1:10, ], type = "response")
#'  intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "probability")
#'  intFit <- predictInterval(gm2, newdata = d1[1:10, ], type = "linear.prediction")
predictInterval <- function(merMod, newdata, level = 0.95,
                            n.sims=100, stat=c("median","mean"),
                            type=c("linear.prediction", "probability"),
                            include.resid.var=TRUE, returnSims = FALSE,
                            seed=NULL, .parallel = FALSE, .paropts = NULL){
  if(any(c("data.frame") != class(newdata))){
    if(any(c("tbl_df", "tbl") %in% class(newdata))){
      newdata <- as.data.frame(newdata)
      warning("newdata is tbl_df or tbl object from dplyr package and has been
              coerced to a data.frame")
    } else{
     newdata <- as.data.frame(newdata)
    }
  }
  predict.type <- match.arg(type,
                            c("linear.prediction", "probability"),
                            several.ok = FALSE)
  stat.type <- match.arg(stat,
                         c("median","mean"),
                         several.ok = FALSE)

  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

  ##First: check if it is a GLMM or NLMM and draw from sigma distribution or incorporate scale parameter if GLMM
  merMod.devcomp <- getME(merMod, "devcomp")
  if (merMod.devcomp$dims[["GLMM"]] == 0 &
      merMod.devcomp$dims[["NLMM"]] == 0) {
    sigmahat <-  sqrt(1/rgamma(n.sims, 0.5 * residDF.merMod(merMod),
                               0.5 * merMod.devcomp$cmp[["pwrss"]]))
    if (predict.type=="probability") {
      predict.type="linear.prediction"
      warning("    Asking for predictions on the probability scale makes no sense, resetting predict.type to linear.prediction",
              call.=FALSE)
    }
  }
  else if (merMod.devcomp$dims[["GLMM"]] == TRUE &
           merMod@resp$family$family == "binomial" &
           merMod@resp$family$link %in% c("logit", "probit")) {
    sigmahat <- rep(1,n.sims)
  }
  else {
    warning("   Prediction for NLMMs or GLMMs that are not mixed binomial regressions is not tested. Sigma set at 1.")
    sigmahat <- rep(1,n.sims)
  }

  matrixForm <- formulaBuild(merMod)
  if(identical(newdata, merMod@frame)){
    newdata.modelMatrix <- model.matrix(matrixForm,
                                        data = merMod@frame)
  } else{
    # combine modelframe trimmed and newdata to ensure proper factor expansion
    tmp <- plyr::rbind.fill(newdata, trimModelFrame(merMod@frame))
    modDV <- as.character(formula(merMod)[2])
    tmp[, modDV] <- 1 # avoid dropping cases without a valid value for the DV
    newdata.modelMatrix <- model.matrix(matrixForm,
                                        data = tmp)[1:nrow(newdata), , drop=FALSE]
    rm(tmp)
  }
  ##Right now I am not multiplying the BLUP variance covariance matrices by our
  ##draw of sigma (for linear models) because their variation is unique.  If anything,
  ##this is where one would multiply them by draws of theta from the model.
  re.xb <- vector(getME(merMod, "n_rfacs"), mode = "list")
  names(re.xb) <- names(ngrps(merMod))
    for(j in names(re.xb)){
    reMeans <- as.matrix(ranef(merMod)[[j]])
    reMatrix <- attr(ranef(merMod, condVar=TRUE)[[j]], which = "postVar")
    # OK, let's knock out all the random effects we don't need
    obslvl <- unique(as.character(newdata[, j]))
    alllvl <- rownames(reMeans)
    keep <- intersect(obslvl, alllvl)
    # Add switch if no random groups are observed to avoid indexing errors,
    # we burn 1 sample of 1 group of all coefficients that will eventually
    # be multiplied by zero later on
    if(length(keep) > 0 & !identical(keep, alllvl)){
      reMeans <- reMeans[keep, , drop=FALSE]
      dimnames(reMatrix)[[3]] <- alllvl
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else if(length(keep) > 0 & identical(keep, alllvl)){
      dimnames(reMatrix)[[3]] <- alllvl
    } else{
      reMeans <- reMeans[1, , drop=FALSE]
      reMatrix <- reMatrix[, , 1, drop = FALSE]
    }
    # -- INSERT chunking code here
    reSimA <- array(data = NA, dim = c(nrow(reMeans), ncol(reMeans), n.sims),
                    dimnames = list(attr(reMeans, "dimnames")[[1]],
                                    attr(reMeans, "dimnames")[[2]],
                                    NULL))
    for(k in 1:nrow(reMeans)){
      meanTmp <- as.matrix(reMeans[k, ])
      matrixTmp <- as.matrix(reMatrix[,,k])
      reSimA[k, ,] <- mvtnorm::rmvnorm(n.sims,
                                       mean=meanTmp,
                                       sigma=matrixTmp,
                                       method="chol") #cholesky is fastest
    }
     tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
     keep <- names(tmp)[names(tmp) %in% dimnames(reSimA)[[2]]]
     tmp <- tmp[, c(keep, "var"), drop = FALSE]
     tmp[, "var"] <- as.character(tmp[, "var"])
     colnames(tmp)[which(names(tmp) == "var")] <- names(newdata[, j,  drop=FALSE])
     tmp.pred <- function(data, coefs, group){
      new.levels <- unique(as.character(data[, group])[!as.character(data[, group]) %in% dimnames(coefs)[[1]]])
       msg <- paste("     The following levels of ", group, " from newdata \n -- ", paste0(new.levels, collapse=", "),
                    " -- are not in the model data. \n     Currently, predictions for these values are based only on the \n fixed coefficients and the observation-level error.", sep="")
       if(length(new.levels > 0)){
         warning(msg, call.=FALSE)
       }
       yhatTmp <- array(data = NA, dim = c(nrow(data), dim(coefs)[3]))
       colIdx <- ncol(data) - 1
       for(i in 1:nrow(data)){
         lvl <- as.character(data[, group][i])
         if(lvl %in% dimnames(coefs)[[1]]){
           yhatTmp[i, ] <- as.numeric(data[i, 1:colIdx]) %*% coefs[lvl, 1:colIdx, ]
         } else{
           # 0 out the RE for these new levels
           yhatTmp[i, ] <- rep(0, colIdx) %*% coefs[1, 1:colIdx, ]
         }
       }
       # for(k in 1:ncol(data)-1){
       #   for(i in 1:nrow(data)){
       #     lvl <- as.character(data[, group][i])
       #     if(lvl %in% dimnames(coefs)[[1]]){
       #       yhatTmp[i, k,] <- data[i, k] * coefs[lvl, k, ]
       #     } else{
       #       yhatTmp[i, k,] <- as.numeric(data[i, k]) * 0
       #     }
       #
       #   }
       # }
       return(yhatTmp)
     }
     # -- INSERT CHUNK COMBINING HERE
     if(nrow(tmp) > 200 | .parallel){
       if(.parallel){
         setup_parallel()
       }
       tmp2 <- split(tmp, (0:nrow(tmp) %/% 100)) #TODO: Find optimum splitting factor
       i <- seq_len(length(tmp2))
       i <- 1:10 # TODO: fix this
       fe_call <- as.call(c(list(quote(foreach::foreach), i = i, .combine = 'rbind', .paropts)))
       fe <- eval(fe_call)
       re.xb[[j]] <- foreach::`%dopar%`(fe, tmp.pred(data = tmp2[[i]],
                                                 coefs =reSimA[, keep, , drop = FALSE],
                                                 group = names(newdata[, j, drop=FALSE])))
     } else{
       re.xb[[j]] <- tmp.pred(data = tmp, coefs = reSimA[, keep, , drop = FALSE],
                              group = names(newdata[, j, drop=FALSE]))
     }

    }
  rm(reSimA)
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
  fe.tmp <- fixef(merMod)
  vcov.tmp <- as.matrix(vcov(merMod))
  if(n.sims > 2000 | .parallel){
    if(.parallel){
      setup_parallel()
    }
    i <- 1:n.sims
    fe_call <- as.call(c(list(quote(foreach::foreach), i = i,
                              .packages = "mvtnorm",
                              .combine = 'rbind', .paropts)))
    fe <- eval(fe_call)
    betaSim <- foreach::`%dopar%`(foreach::foreach(i = 1:n.sims, .combine = "rbind"),
                                  mvtnorm::rmvnorm(1, mean = fe.tmp, sigma = sigmahat[[i]]*vcov.tmp,
                                                   method="chol"))
  } else {
    betaSim <- abind::abind(lapply(1:n.sims,
                                   function(x) mvtnorm::rmvnorm(1, mean = fe.tmp, sigma = sigmahat[x]*vcov.tmp,
                                                                method="chol")), along=1)
  }
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
  rm(re.xb)
  N <- nrow(newdata)
  outs <- data.frame("fit" = rep(NA, N),
                     "upr" = rep(NA, N),
                     "lwr" = rep(NA, N))
  upCI <- 1 - ((1-level)/2)
  loCI <- ((1-level)/2)
  if (include.resid.var==TRUE)
    yhat <- abind::abind(lapply(1:n.sims, function(x) rnorm(N, yhat[,x], sigmahat[x])), along = 2)
  #Output prediction intervals
  if (stat.type == "median") {
    outs[, 1:3] <- t(apply(yhat, 1, quantile, prob = c(0.5, upCI, loCI), na.rm=TRUE))
  }
  if (stat.type == "mean") {
    outs$fit <- apply(yhat, 1, mean, na.rm=TRUE)
    outs[, 2:3] <- t(apply(yhat, 1, quantile, prob = c(upCI, loCI), na.rm=TRUE))
  }
  if (predict.type == "probability") {
    outs <- apply(outs, 2, merMod@resp$family$linkinv)
  }
  #Close it out
  if(returnSims == FALSE){
    return(as.data.frame(outs))
  } else if(returnSims == TRUE){
    outs <- as.data.frame(outs)
    attr(outs, "sim.results") <- yhat
    return(outs)
  }
}
