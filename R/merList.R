# Imputation functions
#' Extract model information from a merMod
#'
#' @param object a merMod object
#'
#' @return Simple summary information about the object, coefficients, number
#' of observations, and adjusted R squared
modelInfo <- function(object){
  out <- rbind(summary(object)$coefficients,
               round(length(object$fitted.values), digits = 0),
               round(summary(object)$adj.r.squared, digits = 5))
  row.names(out) <- c(row.names(summary(object)$coefficients), "n.obs", "adj. rsquared")
  return(out)
}

# Functions to extract standard deviation of random effects from model
#' Extract the standard deviation of the random effects from a merMod object
#'
#' @param model an object that inherits from class merMod
#'
#' @return a numeric vector for standard deviations of the random effects
#' @export
REsdExtract <- function(model){
  out <- unlist(lapply(VarCorr(model), attr, "stddev"))
  return(out)
}

#' Extract the correlations between the slopes and the intercepts from a model
#'
#' @param model an object that inherits from class merMod
#'
#' @return a numeric vector of the correlations among the effects
#' @export
REcorrExtract <- function(model){
  out <- unlist(lapply(VarCorr(model), attr, "corre"))
  return(min(unique(out)))
}


#' Extract data.frame of random effect statistics from merMod List
#'
#' @param modList a list of multilevel models
#'
#' @return a data.frame
#' @import dplyr
#' @export
modelRandEffStats <- function(modList){
  effList <- lapply(modList, tidy, effects = "ran_pars")
  effList <- do.call(rbind, effList)
  out <- effList %>% group_by(term, group) %>%
    summarize(est = mean(estimate),
              std.error = sd(estimate)) %>%
    rename(estimate = est)
  return(as.data.frame(out))
}

#' Extract averaged fixed effect parameters across a list of merMod objects
#'
#' @param modList an object of class merModList
#' @param ... additional arguments to pass to \code{\link{tidy}}
#'
#' @return a data.frame of the averaged fixed effect parameters
#' @export
#' @importFrom broom tidy
#' @import dplyr
modelFixedEff <- function(modList, ...){
  fixEst <- lapply(modList, tidy, effects = "fixed", ...)
  fixEst <- do.call(rbind, fixEst)
  # Collapse
  out <- fixEst %>% dplyr::group_by(term) %>%
          dplyr::summarize(estimate = mean(estimate),
               std.error = mean(std.error))
  out$statistic <- out$estimate / out$std.error
  return(as.data.frame(out))
}

#' Extract the variances and correlations for random effects from a merMod list
#' @inheritParams lme4::VarCorr
#' @param rdig the number of digits to round to, integer
#' @return a list with two elements "stddev" and "correlation" for the standard
#' deviations and correlations averaged across models in the list
#' @export
#' @import lme4
VarCorr.merModList <- function(x, sigma = 1, rdig = 3L){
  modList <- x
  ngrps <- length(VarCorr(modList[[1]]))
  errorList <- vector(mode = 'list', length = ngrps)
  corrList <- vector(mode = 'list', length = ngrps)
  for(i in 1:ngrps){
    subList <- lapply(modList, function(x) VarCorr(x)[[i]])
    if(all(dim(subList[[1]]) == c(1, 1))){
      subList <- mean(sqrt(unlist(subList)))
      errorList[[i]] <- subList
      names(errorList) <- "Intercept"
      corrList[[i]] <- matrix(1)
      dimnames(corrList[[i]]) <- list("(Intercept)", "(Intercept)")
    } else {
      errorList[[i]] <- apply(simplify2array(lapply(subList, attr, "stddev")),
                              1, mean)
      corrList[[i]] <- apply(simplify2array(lapply(subList, attr, "corre")),
                             1:2,mean)
    }
  }
  for(i in 1:length(errorList)){
    if(is.null(names(errorList[[i]]))){
      names(errorList[[i]]) <- "(Intercept)"
    }
  }
  for(i in 1:length(corrList)){
    if(is.null(names(corrList[[i]])) & is.null(dim(corrList[[i]]))){
      names(corrList[[i]]) <- "(Intercept)"
    }
  }
  names(errorList) <- names(ranef(modList[[1]]))
  names(corrList) <- names(ranef(modList[[1]]))
  return(list("stddev" = errorList, "correlation" = corrList))
}

utils::globalVariables(c("term", "estimate","std.error"))
#' Print the results of a merMod list
#'
#' @param x a modList of class merModList
#' @param ... additional arguments
#'
#' @return summary content printed to console
#' @export
print.merModList <- function(x, ...){
  modList <- x
  args <- eval(substitute(alist(...)))
  if("digits" %in% names(args)){
    digits <- args$digits
  } else{
    digits <- 3
  }
  len <- length(modList)
  form <- modList[[1]]@call
  print(summary(modList[[1]])$methTitle)
  cat("Model family: ", summary(modList[[1]])$family)
  cat("\n")
  print(form)
  cat("\nFixed Effects:\n")
  fedat <- modelFixedEff(modList)
  dimnames(fedat)[[1]] <- fedat$term
  pfround(fedat[, -1], digits)
  cat("\nRandom Effects:\n")
  ngrps <- length(VarCorr(modList[[1]]))
  errorList <- VarCorr(modList)$stddev
  corrList <- VarCorr(modList)$correlation

  cat("\nError Term Standard Deviations by Level:\n")
  for(i in 1:length(errorList)){
    cat("\n")
    cat(names(errorList[i]))
    cat("\n")
    if(is.null(names(errorList[[i]]))){
      names(errorList[[i]]) <- "(Intercept)"
    }
    pfround(errorList[[i]], digits = digits)
    cat("\n")
  }
  # lapply(errorList, pfround, digits)
  cat("\nError Term Correlations:\n")
  for(i in 1:length(corrList)){
    cat("\n")
    cat(names(corrList[i]))
    cat("\n")
    if(is.null(names(corrList[[i]]))){
      names(corrList[[i]]) <- "(Intercept)"
    }
    pfround(corrList[[i]], digits = digits)
    cat("\n")
  }
  # lapply(corrList, pfround, digits)
  residError <- mean(unlist(lapply(modList, function(x) attr(VarCorr(x), "sc"))))
  cat("\nResidual Error =", fround(residError,
                                   digits), "\n")
  cat("\n---Groups\n")
  ngrps <- lapply(modList[[1]]@flist, function(x) length(levels(x)))
  modn <- getME(modList[[1]], "devcomp")$dims["n"]
  cat(sprintf("number of obs: %d, groups: ", modn))
  cat(paste(paste(names(ngrps), ngrps, sep = ", "),
            collapse = "; "))
  cat("\n")
  cat("\nModel Fit Stats")
  mAIC <- mean(unlist(lapply(modList, AIC)))
  cat(sprintf("\nAIC = %g", round(mAIC, 1)))
  moDsigma.hat <- mean(unlist(lapply(modList, sigma)))
  cat("\nResidual standard deviation =", fround(moDsigma.hat,
                                                digits), "\n")
}


#' Apply a multilevel model to a list of data frames
#'
#' @param formula a formula to pass through compatible with merMod
#' @param data a list object with each element being a data.frame
#' @param parallel logical, should the models be run in parallel?
#' @param ... additional arguments to pass to the estimating function
#' @rdname merModList
#'
#' @return a list of fitted merMod objects of class merModList
#' @export
lmerModList <- function(formula, data, parallel = NULL, ...){
  ml <- lapply(data, function(d) lmer(formula, data = d, ...))
  class(ml) <- "merModList"
  return(ml)
}

#' Apply a Bayesian multilevel model to a list of data frames
#'
#' @inheritParams lmerModList
#' @rdname merModList
#' @return a merModList
#' @importFrom blme blmer
#' @export
blmerModList <- function(formula, data, parallel = NULL, ...){
  ml <- lapply(data, function(d) blmer(formula, data = d, ...))
  class(ml) <- "merModList"
  return(ml)
}

#' Apply a generalized linear multilevel model to a list of data frames
#'
#' @inheritParams lmerModList
#' @rdname merModList
#' @return a merModList
#' @export
glmerModList <- function(formula, data, parallel = NULL, ...){
  ml <- lapply(data, function(d) glmer(formula, data = d, ...))
  class(ml) <- "merModList"
  return(ml)
}

#' Apply a Bayesian generalized linear multilevel model to a list of data frames
#'
#' @inheritParams lmerModList
#' @rdname merModList
#' @return a merModList
#' @importFrom blme bglmer
#' @export
bglmerModList <- function(formula, data, parallel = NULL, ...){
  ml <- lapply(data, function(d) bglmer(formula, data = d, ...))
  class(ml) <- "merModList"
  return(ml)
}

