#' Predict from merMod objects with a prediction interval
#' @description This function provides a way to capture model uncertainty in
#' predictions from multi-level models fit with \code{lme4}. By drawing a sampling
#' distribution for the random and the fixed effects and then estimating the fitted
#' value across that distribution, it is possible to generate a prediction interval
#' for fitted values that includes all variation in the model except for variation
#' in the covariance parameters, theta. This is a much faster alternative than
#' bootstrapping for models fit to medium to large datasets.
#' @param merMod a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param which a character specifying what to return, by default it returns the
#' full interval, but you can also select to return only the fixed variation or
#' the random component variation. If full is selected the resulting data.frame
#' will be \code{nrow(newdata) * number of model levels} long
#' @param level the width of the prediction interval
#' @param n.sims number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param type type of prediction to develop
#' @param include.resid.var logical, include or exclude the residual variance for
#' linear models
#' @param returnSims logical, should all n.sims simulations be returned?
#' @param seed numeric, optional argument to set seed for simulations
#' @param fix.intercept.variance logical; should the variance of the intercept
#' term be adjusted downwards to roughly correct for its covariance with the
#' random effects, as if all the random effects are intercept effects?
#' @param ignore.fixed.terms a numeric or string vector of indexes or names of
#' fixed effects which should be considered as fully known (zero variance). This
#' can result in under-conservative intervals, but for models with random effects
#' nested inside fixed effects, holding the fixed effects constant intervals may
#' give intervals with closer to nominal coverage than the over-conservative
#' intervals without this option, which ignore negative correlation between the
#' outer (fixed) and inner (random) coefficients.
#' @param .parallel, logical should parallel computation be used, default is FALSE
#' @param .paropts, -NOT USED: Caused issue #54- a list of additional options passed into the foreach function
#' when parallel computation is enabled. This is important if (for example) your
#' code relies on external data or packages: use the .export and .packages arguments
#'  to supply them so that all cluster nodes have the correct environment set up
#'  for computing.
#' @return a data.frame with three columns:
#' \describe{
#'     \item{\code{fit}}{The center of the distribution of predicted values as defined by
#'     the \code{stat} parameter.}
#'     \item{\code{lwr}}{The lower prediction interval bound corresponding to the quantile cut
#'     defined in \code{level}.}
#'     \item{\code{upr}}{The upper prediction interval bound corresponding to the quantile cut
#'     defined in \code{level}.}
#'   }
#' If returnSims = TRUE, then the individual simulations are attached to this
#' data.frame in the attribute \code{sim.results} and are stored as a matrix.
#' @details To generate a prediction interval, the function first computes a simulated
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
#' @import lme4
#' @importFrom abind abind
#' @importFrom mvtnorm rmvnorm
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @examples
#' \donttest{
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
#'  }
predictInterval <- function(merMod, newdata, which=c("full", "fixed", "random", "all"),
                            level = 0.8,
                            n.sims = 1000, stat=c("median","mean"),
                            type=c("linear.prediction", "probability"),
                            include.resid.var=TRUE, returnSims = FALSE,
                            seed=NULL, .parallel = FALSE, .paropts = NULL,
                            fix.intercept.variance = FALSE, #This does NOT work with random slope models
                            ignore.fixed.terms = NULL)
                            {
  if(missing(newdata)){
    newdata <- merMod@frame
  }
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
  which.eff <- match.arg(which,
                         c("full", "fixed", "random", "all"),
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

  newdata.modelMatrix <- buildModelMatrix(model= merMod, newdata = newdata)
  # When there is no fixed effect intercept but there is a group level intercept
  # We need to do something!

  rr <- ranef(merMod, condVar = TRUE)
  re.xb <- vector(getME(merMod, "n_rfacs"), mode = "list")
  names(re.xb) <- names(ngrps(merMod))
  for (j in names(re.xb)){
    reMeans <- as.matrix(rr[[j]])
    reMatrix <- attr(rr[[j]], which = "postVar")
    # OK, let's knock out all the random effects we don't need
    if (j %in% names(newdata)){ # get around if names do not line up because of nesting
      obslvl <- unique(as.character(newdata[, j]))
      alllvl <- rownames(reMeans)
      keep <- intersect(obslvl, alllvl)
    } else {
      obslvl <- colnames(newdata.modelMatrix)
      alllvl <- rownames(reMeans)
      keep <- intersect(obslvl, alllvl)
    }
    # Add switch if no random groups are observed to avoid indexing errors,
    # we burn 1 sample of 1 group of all coefficients that will eventually
    # be multiplied by zero later on
    if (length(keep) > 0 & !identical(keep, alllvl)) {
      reMeans <- reMeans[keep, , drop=FALSE]
      dimnames(reMatrix)[[3]] <- alllvl
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else if (length(keep) > 0 & identical(keep, alllvl)){
      dimnames(reMatrix)[[3]] <- alllvl
      # dimnames(reMeans)[[2]] <- j  # we need to get the variable name into this ojbect
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else{
      reMeans <- reMeans[1, , drop=FALSE]
      reMatrix <- reMatrix[, , 1, drop = FALSE]
    }

    tmpList <- vector(length = nrow(reMeans), mode = "list")
    for (k in 1:nrow(reMeans)){
      meanTmp <- reMeans[k, ]
      names(meanTmp) <- NULL
      matrixTmp <- as.matrix(reMatrix[, , k])
      tmpList[[k]] <- as.matrix(mvtnorm::rmvnorm(n= n.sims,
                                                mean=meanTmp,
                                                sigma=matrixTmp, method = "chol"))
    }

    REcoefs <- sapply(tmpList, identity, simplify="array")
    # rm(tmpList)
    dimnames(REcoefs) <- list(1:n.sims,
                            attr(reMeans, "dimnames")[[2]],
                            attr(reMeans, "dimnames")[[1]]
                            )
    if (j %in% names(newdata)) { # get around if names do not line up because of nesting
      newdata.modelMatrix <- as.matrix(newdata.modelMatrix)  ## give up, sparse to dense now
      tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
      tmp <- tmp[, !duplicated(colnames(tmp))]
      keep <- names(tmp)[names(tmp) %in% colnames(REcoefs)]
      if (length(keep) == 0) {
        keep <- grep(dimnames(REcoefs)[[2]], names(tmp), value = TRUE)
      }
        if (length(keep) == 0) {
          tmp <- cbind(model.frame(subbars(formula(merMod)), data = newdata),
                       var = newdata[, j])
          keep <- grep(dimnames(REcoefs)[[2]], names(tmp), value = TRUE)
        }
          if ( length(keep) == 0) {
            # Add in an intercept for RE purposes
            tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
            tmp <- tmp[, !duplicated(colnames(tmp))]
            tmp <- cbind(data.frame(1), tmp)
            names(tmp)[1] <- "(Intercept)"
            keep <- "(Intercept)"
          }
      tmp <- tmp[, c(keep, "var"), drop = FALSE]
      tmp[, "var"] <- as.character(tmp[, "var"])
      colnames(tmp)[which(names(tmp) == "var")] <- names(newdata[, j,  drop = FALSE])
      if (all(grepl(":", keep))) {
        # Strip out the interaction after
        keep <- unique(gsub("(.*):.*", "\\1", keep))
      }
    } else {
      tmp <- as.data.frame(newdata.modelMatrix)
      tmp <- tmp[, !duplicated(colnames(tmp))] # deduplicate columns because
      # column names can be duplicated to account for multiple effects
      # but we've already reconciled all the effects
      tmp$var <- names(tmp[keep])[max.col(tmp[keep])] #changed alllvl to keep in
      #this line re: issue #53 where newdata doesn't have all levels of rfx in
      #nested specification (with ":") so this just takes the subset of alllvl
      #that are specified in model
      keep <- names(tmp)[names(tmp) %in% dimnames(REcoefs)[[2]]]
      tmp <- tmp[, c(keep, "var"), drop = FALSE]
      tmp[, "var"] <- as.character(tmp[, "var"])
      colnames(tmp)[which(names(tmp) == "var")] <- j
    }
    #######################
    ################
     tmp.pred <- function(data, coefs, group){
      new.levels <- unique(as.character(data[, group])[!as.character(data[, group]) %in% dimnames(coefs)[[3]]])
       msg <- paste("     The following levels of ", group, " from newdata \n -- ", paste0(new.levels, collapse=", "),
                    " -- are not in the model data. \n     Currently, predictions for these values are based only on the \n fixed coefficients and the observation-level error.", sep="")
       if(length(new.levels > 0)){
         warning(msg, call.=FALSE)
       }
       yhatTmp <- array(data = NA, dim = c(nrow(data), dim(coefs)[1]))
       colIdx <- ncol(data) - 1

      colLL <- length(1:colIdx)
      if(colLL > dim(coefs)[2]) {
        # copy over
        coefs_new <- array(NA, dim = c(dim(coefs)[1], colLL,
                                       dim(coefs)[3]))
        dimnames(coefs_new)[c(1, 3)] <- dimnames(coefs)[c(1, 3)]

        dimnames(coefs_new)[[2]] <- rep(dimnames(coefs)[[2]], dim(coefs_new)[2])
        for (k in 1:colLL) {
          coefs_new[, k, 1:dim(coefs)[3]] <- coefs[, 1, 1:dim(coefs)[3]]
        }
        coefs <- coefs_new
      }

      for(i in 1:nrow(data)){
        lvl <- as.character(data[, group][i])
         if(!lvl %in% new.levels){
           yhatTmp[i, ] <- as.numeric(data[i, 1:colIdx]) %*% t(coefs[, 1:colIdx, lvl])
         } else{
           # 0 out the RE for these new levels
           yhatTmp[i, ] <- rep(0, colIdx) %*% t(coefs[, 1:colIdx, 1])
         }
       }
       rownames(yhatTmp) <- rownames(data)
       rm(data)
       return(yhatTmp)
     }
    #########################
    ####
     if(nrow(tmp) > 1000 | .parallel) {
       if (requireNamespace("foreach", quietly=TRUE)) {
         if(.parallel){
         setup_parallel()
       }
       tmp2 <- split(tmp, (1:nrow(tmp) %/% 500)) #TODO: Find optimum splitting factor
       tmp2 <- tmp2[lapply(tmp2,length) > 0]
       fe_call <- as.call(c(list(quote(foreach::foreach), i = seq_along(tmp2),
                                 .combine = 'rbind')))
       fe <- eval(fe_call)
       re.xb[[j]] <- foreach::`%dopar%`(fe, tmp.pred(data = tmp2[[i]],
                                                 coefs = REcoefs[, keep, , drop = FALSE],
                                                 group = j))
       rm(tmp2)
       } else {
         warning("foreach package is unavailable, parallel computing not available")
         re.xb[[j]] <- tmp.pred(data = tmp, coefs = REcoefs[, keep, , drop = FALSE],
                                group = j)
       }
     } else{
       re.xb[[j]] <- tmp.pred(data = tmp, coefs = REcoefs[, keep, , drop = FALSE],
                              group = j)
     }
     rm(tmp)
    }
  rm(REcoefs)
  # TODO: Add a check for new.levels that is outside of the above loop
  # for now, ignore this check
  if (include.resid.var==FALSE) {
#     if (length(new.levels)==0)
      sigmahat <- rep(1, n.sims)
#     else {
#       include.resid.var=TRUE
#       warning("    \n  Since new levels were detected resetting include.resid.var to TRUE.")
#     }
  }

  # fixed.xb is nrow(newdata) x n.sims
  ##Calculate yhat as sum of the components (fixed plus all groupling factors)
  fe.tmp <- fixef(merMod)
  vcov.tmp <- as.matrix(vcov(merMod))

  # Detect if an intercept is present
  # TODO - is this reliable
  if (is.na(names(attr(VarCorr(merMod)[[j]],"stddev")["(Intercept)"]))) {
    fix.intercept.variance <- FALSE
    message("No intercept detected, setting fix.intercept.variance to FALSE")
  }
  # If intercept is not in fixed terms
  if (!"(Intercept)" %in% names(fixef(merMod)) && fix.intercept.variance) {
    # TODO - decide if this is an error or if we should allow it to continue with warning
    warning("No fixed-effect intercept detected. Variance adjustment may be unreliable.")
  }

  if (fix.intercept.variance) {
    #Assuming all random effects include intercepts.
    intercept.variance <- vcov.tmp[1,1]

    groupsizes <- ngrps(merMod)
    for(j in names(groupsizes)){ #for every group of random e

      groupExtraPrecision <- 0
      groupVar <- (attr(VarCorr(merMod)[[j]],"stddev")["(Intercept)"])^2
      reMatrix <- attr(rr[[j]], which = "postVar")
      for (eff in 1:dim(reMatrix)[3]) {
        term <- 1/(reMatrix[1,1,eff] + groupVar)
        if (term > 0) {
            groupExtraPrecision <- groupExtraPrecision + term
        } else {
            warning("fix.intercept.variance got negative precision; better turn it off.")
        }
      }
      intercept.variance <- intercept.variance - 1/groupExtraPrecision
    }


    if (intercept.variance < 0) {
        warning("fix.intercept.variance got negative variance; better turn it off.")
    }
    ratio <- intercept.variance/vcov.tmp[1,1]
    prec.tmp <- solve(vcov.tmp)
    prec.tmp[1,1] <- prec.tmp[1,1] / ratio
    vcov.tmp[1,] <- vcov.tmp[1,] * ratio
    vcov.tmp <- solve(prec.tmp, tol=1e-50)
  }
  if (!is.null(ignore.fixed.terms)) {
      prec.tmp <- solve(vcov.tmp)
      for (term in ignore.fixed.terms) {
            prec.tmp[term,term] <- prec.tmp[term,term] * 1e15
      }
      vcov.tmp <- solve(prec.tmp, tol=1e-50)
  }
  if(n.sims > 2000 | .parallel){
    if(.parallel){
      setup_parallel()
    }
    i <- 1:n.sims
    fe_call <- as.call(c(list(quote(foreach::foreach), i = i,
                              .combine = 'rbind')))
    fe <- eval(fe_call)
    betaSim <- foreach::`%dopar%`(fe, mvtnorm::rmvnorm(n = 1, mean = fe.tmp,
                                                sigma = vcov.tmp,
                                  method = "chol"))

  } else {
    betaSim <- abind::abind(lapply(1:n.sims,
                               function(x) mvtnorm::rmvnorm(n = 1, mean = fe.tmp,
                                                    sigma = vcov.tmp,
                                method = "chol")), along=1)
  }
  # Pad betaSim
  colnames(betaSim) <- names(fe.tmp)
  rownames(betaSim) <- 1:n.sims
  newdata.modelMatrix <- buildModelMatrix(merMod, newdata = newdata, which = "fixed")
  if (ncol(newdata.modelMatrix) > ncol(betaSim)) {
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
  ######
  if(which.eff == "full"){
    yhat <- Reduce('+', re.xb)
  } else if(which.eff == "fixed"){
    yhat <- Reduce('+', re.xb["fixed"])
  } else if(which.eff == "random"){
    re.xb["fixed"] <- NULL
    yhat <- Reduce('+', re.xb)
  } else if(which.eff == "all"){
    yhat <- Reduce('+', re.xb)
    N <- nrow(newdata)
    if (include.resid.var==TRUE){
      for(i in 1:length(re.xb)){
        re.xb[[i]] <- abind::abind(lapply(1:n.sims, function(x) rnorm(N, re.xb[[i]][, x], sigmahat[x])), along=2)
      }
    }
    pi.comps <- re.xb
  }
  rm(re.xb)
  N <- nrow(newdata)
  outs <- data.frame("fit" = rep(NA, N),
                     "upr" = rep(NA, N),
                     "lwr" = rep(NA, N))
  upCI <- 1 - ((1-level)/2)
  loCI <- ((1-level)/2)
  if (include.resid.var==TRUE){
    yhat <- abind::abind(lapply(1:n.sims,
                                function(x) rnorm(N, yhat[,x], sigmahat[x])),
                         along = 2)
  }
  # Output prediction intervals
  if (stat.type == "median") {
    outs[, 1:3] <- t(apply(yhat, 1, quantile, prob = c(0.5, upCI, loCI),
                           na.rm=TRUE))
  }
  if (stat.type == "mean") {
    outs$fit <- apply(yhat, 1, mean, na.rm=TRUE)
    outs[, 2:3] <- t(apply(yhat, 1, quantile, prob = c(upCI, loCI),
                           na.rm=TRUE))
  }
  if (predict.type == "probability") {
    if(nrow(outs) == 1) {
      outs <- t(apply(outs, 2, merMod@resp$family$linkinv))
    } else {
      outs <- apply(outs, 2, merMod@resp$family$linkinv)
    }

  }

  ##############################
  # Construct observation predictors for each component of the model
  ##########################
  if(which.eff == "all"){
    if(returnSims == TRUE){
      allSims <- pi.comps
    }
    for(i in 1:length(pi.comps)){
      if( stat.type == "median"){
        pi.comps[[i]] <-  t(apply(pi.comps[[i]], 1, quantile,
                                  prob = c(0.5, upCI, loCI), na.rm=TRUE))
        pi.comps[[i]] <- as.data.frame(pi.comps[[i]])
        names(pi.comps[[i]]) <- c("fit", "upr", "lwr")
      }
      if(stat.type == "mean"){
        tmp <- pi.comps[[i]]
        pi.comps[[i]] <- data.frame("fit" = rep(NA, N), "upr" =NA,
                          "lwr" = NA)
        pi.comps[[i]]$fit <- apply(tmp, 1, mean, na.rm=TRUE)
        pi.comps[[i]][, 2:3] <- t(apply(tmp, 1, quantile, prob = c(upCI, loCI), na.rm=TRUE))
      }
      if (predict.type == "probability") {
        pi.comps[[i]] <- apply(pi.comps[[i]], 2, merMod@resp$family$linkinv)
        pi.comps[[i]] <- as.data.frame(pi.comps[[i]])
        names(pi.comps[[i]]) <- c("fit", "upr", "lwr")
      }
    }
    componentOut <- dplyr::bind_rows(pi.comps, .id="effect")
    outs <- cbind(data.frame("effect" = "combined"), outs)
    outs <- suppressWarnings(bind_rows(outs, componentOut))
    outs$obs <- rep(1:N, nrow(outs) %/% N)
    rm(pi.comps)
  }

  #Close it out
  if(returnSims == FALSE){
    return(as.data.frame(outs))
  } else if(returnSims == TRUE){
    outs <- as.data.frame(outs)
    if(which.eff == "all"){
      attr(outs, "sim.results") <- allSims
    } else{
      attr(outs, "sim.results") <- yhat
    }
    return(outs)
  }
}

## TODO: Finish exporting so that all returns the individual predictions for
# each random effect separately
