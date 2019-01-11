# Imputation functions
#' Extract model information from a merMod
#'
#' @param object a merMod object
#'
#' @return Simple summary information about the object, number
#' of observations, number of grouping terms, AIC, and residual standard deviation
#' @export
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' modelInfo(mod[[1]])
#' lapply(mod, modelInfo)
modelInfo <- function(object){
  ngrps <- lapply(object@flist, function(x) length(levels(x)))
  out <- data.frame("n.obs" = getME(object, "devcomp")$dims["n"],
                    "n.lvls" = length(ngrps),
                    "AIC" = AIC(object),
                    "sigma" = sigma(object))
  row.names(out) <- NULL
  # cat("\n---Groups\n")
  # ngrps <- lapply(modList[[1]]@flist, function(x) length(levels(x)))
  # modn <- getME(modList[[1]], "devcomp")$dims["n"]
  # cat(sprintf("number of obs: %d, groups: ", modn))
  # cat(paste(paste(names(ngrps), ngrps, sep = ", "),
  #           collapse = "; "))
  # cat("\n")
  # cat("\nModel Fit Stats")
  # mAIC <- mean(unlist(lapply(modList, AIC)))
  # cat(sprintf("\nAIC = %g", round(mAIC, 1)))
  # moDsigma.hat <- mean(unlist(lapply(modList, sigma)))
  # cat("\nResidual standard deviation =", fround(moDsigma.hat,
  #                                               digits), "\n")
  return(out)

}

# Functions to extract standard deviation of random effects from model
#' Extract the standard deviation of the random effects from a merMod object
#'
#' @param model an object that inherits from class merMod
#'
#' @return a numeric vector for standard deviations of the random effects
#' @export
#' @examples
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' REsdExtract(fm1)
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
#' @examples
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' REcorrExtract(fm1)
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
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' modelRandEffStats(mod)
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
#' @details The Rubin correction for combining estimates and standard errors from
#' Rubin (1987) is applied to adjust for the within and between imputation variances.
#' @export
#' @importFrom broom tidy
#' @import dplyr
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' modelFixedEff(mod)
modelFixedEff <- function(modList, ...){
  fixEst <- lapply(modList, tidy, effects = "fixed", ...)
  fixEst <- do.call(rbind, fixEst)
  # Collapse
  # Rubin correction, get length of list
  ml <- length(modList)
  # Get between and within imputation variance, apply total correction
  # Calculate degree of freedom correction
  rubin <- fixEst %>% group_by(term) %>%
    mutate(mean_est = mean(estimate)) %>%
    mutate(est_ss = (estimate - mean_est)^2) %>%
    summarize(estimate = mean(estimate),
              within_var = mean(std.error), # compute within imputation variance
              between_var =  mean(est_ss)) %>% # estimate the between imputation variance
    mutate(std.error = within_var + ((1 + 1/ml)*between_var),
           df = (ml-1)* (1 + within_var/((1 + 1/ml)*between_var))^2) # apply rubins total variance correction
  # add fallback
  if (any((((1 + 1/ml)*rubin$between_var)^2) < 0.000000001)) {
    warning("Between imputation variance is very small, are imputation sets too similar?")
  }
  # DEPRECATED method
  # out <- fixEst %>% dplyr::group_by(term) %>%
  #         dplyr::summarize(estimate = mean(estimate),
  #              std.error = mean(std.error))
  rubin$statistic <- rubin$estimate / rubin$std.error
  rubin <- rubin %>% dplyr::select(term, estimate, std.error, statistic, df)
  return(as.data.frame(rubin))
}


#' Extract fixed-effects estimates for a merModList
#'
#' @inheritParams lme4::fixef
#' @return a named, numeric vector of fixed-effects estimates.
#' @details Extract the estimates of the fixed-effects parameters from a list of
#' fitted \code{merMod} models. Takes the mean of the individual \code{fixef}
#' objects for each of the component models in the \code{merModList}.
#' @export
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' fixef(mod)
fixef.merModList <- function(object, add.dropped = FALSE, ...){
  Reduce(`+`, lapply(object, fixef)) / length(object)
}

#' Extract random-effects estimates for a merModList
#'
#' @inheritParams lme4::ranef
#' @return a named, numeric vector of random-effects estimates.
#' @details Extract the estimates of the random-effects parameters from a list of
#' fitted \code{merMod} models. Takes the mean of the individual \code{ranef}
#' objects for each of the component models in the \code{merModList}.
#' @export
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' ranef(mod)
ranef.merModList <- function(object, ...){
  levels <- getME(object[[1]], "n_rfacs")
  re <- vector(length = levels, mode = "list")
  for(i in seq_along(1:levels)){
    # <- Reduce(`+`, lapply(object, ranef)[i]) / length(object)
  re[i]  <- lapply(Reduce(`+`, lapply(object, ranef)[1]), function(x) x/length(object))
  }
  names(re) <- names(ranef(object[[1]]))
  return(re)
}


#' Extract the variances and correlations for random effects from a merMod list
#' @inheritParams lme4::VarCorr
#' @param rdig the number of digits to round to, integer
#' @return a list with two elements "stddev" and "correlation" for the standard
#' deviations and correlations averaged across models in the list
#' @export
#' @import lme4
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' VarCorr(mod)
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
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' print(mod)
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

#' Summarize a merMod list
#'
#' @param object a modList of class merModList
#' @param ... additional arguments
#'
#' @return a summary object of model information
#' @export
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' summary(mod)
summary.merModList <- function(object, ...){
  out <- lapply(object, sum.mm)
  class(out) <- "summary.merModList"
  return(out)
}

#' Print the summary of a merMod list
#'
#' @param x a summary of amerModList object
#' @param ... additional arguments
#'
#' @return summary content printed to console
#' @export
print.summary.merModList <- function(x, ...){
  lapply(x, print)
}

#' Apply a multilevel model to a list of data frames
#'
#' @param formula a formula to pass through compatible with merMod
#' @param data a list object with each element being a data.frame
#' @param parallel logical, should the models be run in parallel? Default FALSE. If so,
#' the `future_lapply` function from the `future.apply` package is used. See
#' details.
#' @param ... additional arguments to pass to the estimating function
#' @rdname merModList
#'
#' @details Parallel computing is provided by the `futures` package, and its
#' extension the `future.apply` package to provide the `future_lapply` function
#' for easy parallel computations on lists. To use this package, simply register
#' a parallel backend using the `plan()` function from `futures` - an example
#' is to use `plan(multisession)`
#'
#' @return a list of fitted merMod objects of class merModList
#' @export
#' @examples
#' sim_list <- replicate(n = 10,
#'         expr = sleepstudy[sample(row.names(sleepstudy), 180),],
#'         simplify=FALSE)
#' fml <- "Reaction ~ Days + (Days | Subject)"
#' mod <- lmerModList(fml, data = sim_list)
#' summary(mod)
lmerModList <- function(formula, data, parallel = FALSE, ...){
  if(parallel) {
    if (requireNamespace("future.apply", quietly=TRUE)) {
      ml <- future.apply::future_lapply(data, function(d) lmer(formula, data = d, ...))
    }
    warning("Parallel set but future.apply not available. Running sequentially.")
    ml <- lapply(data, function(d) lmer(formula, data = d, ...))
  } else {
    ml <- lapply(data, function(d) lmer(formula, data = d, ...))
  }

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
blmerModList <- function(formula, data, parallel = FALSE, ...){
  if(parallel) {
    if (requireNamespace("future.apply", quietly=TRUE)) {
      ml <- future.apply::future_lapply(data, function(d) blmer(formula, data = d, ...))
    }
    warning("Parallel set but future.apply not available. Running sequentially.")
    ml <- lapply(data, function(d) blmer(formula, data = d, ...))
  } else {
    ml <- lapply(data, function(d) blmer(formula, data = d, ...))
  }
  class(ml) <- "merModList"
  return(ml)
}

#' Apply a generalized linear multilevel model to a list of data frames
#'
#' @inheritParams lmerModList
#' @rdname merModList
#' @return a merModList
#' @export
glmerModList <- function(formula, data, parallel = FALSE, ...){
  if(parallel) {
    if (requireNamespace("future.apply", quietly=TRUE)) {
      ml <- future.apply::future_lapply(data, function(d) glmer(formula, data = d, ...))
    }
    warning("Parallel set but future.apply not available. Running sequentially.")
    ml <- lapply(data, function(d) glmer(formula, data = d, ...))
  } else {
    ml <- lapply(data, function(d) glmer(formula, data = d, ...))
  }
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
bglmerModList <- function(formula, data, parallel = FALSE, ...){
  if(parallel) {
    if (requireNamespace("future.apply", quietly=TRUE)) {
      ml <- future.apply::future_lapply(data, function(d) bglmer(formula, data = d, ...))
    }
    warning("Parallel set but future.apply not available. Running sequentially.")
    ml <- lapply(data, function(d) bglmer(formula, data = d, ...))
  } else {
    ml <- lapply(data, function(d) bglmer(formula, data = d, ...))
  }
  class(ml) <- "merModList"
  return(ml)
}

