#' Predict from merMod objects with a prediction interval
#'
#' @description This function provides a way to capture model uncertainty in
#' predictions from multi-level models fit with \code{lme4}. By drawing a sampling
#' distribution for the random and the fixed effects and then estimating the fitted
#' value across that distribution, it is possible to generate a prediction interval
#' for fitted values that includes all variation in the model except for variation
#' in the covariance parameters, theta. This is a much faster alternative than
#' bootstrapping for models fit to medium to large datasets.
#'
#' @param merMod a merMod object from lme4
#' @param newdata a data.frame of new data to predict
#' @param which a character specifying what to return, by default it returns the
#'   full interval, but you can also select to return only the fixed variation or
#'   the random component variation. If full is selected the resulting data.frame
#'   will be \code{nrow(newdata) * number of model levels} long
#' @param level the width of the prediction interval
#' @param n.sims number of simulation samples to construct
#' @param stat take the median or mean of simulated intervals
#' @param type type of prediction to develop
#' @param include.resid.var logical, include or exclude the residual variance for
#'   linear models
#' @param returnSims logical, should all n.sims simulations be returned?
#' @param seed numeric, optional argument to set seed for simulations
#' @param .parallel logical, use parallel computation (default FALSE)
#' @param .paropts not used - placeholder for future foreach options
#' @param fix.intercept.variance logical; should the variance of the intercept
#'   term be adjusted downwards to roughly correct for its covariance with the
#'   random effects, as if all the random effects are intercept effects?
#' @param ignore.fixed.terms a numeric or string vector of indexes or names of
#'   fixed effects which should be considered as fully known (zero variance).
#'
#' @return a data.frame with three columns: fit, lwr and upr. If `returnSims`
#'   is TRUE the attribute `sim.results` contains the full simulation array.
#' @export
#' @import lme4
#' @importFrom abind abind
#' @importFrom mvtnorm rmvnorm
#' @importFrom foreach %dopar% foreach
predictInterval <- function(
  merMod,
  newdata = NULL,
  which = c("full", "fixed", "random", "all"),
  level = 0.8,
  n.sims = 1000,
  stat = c("median", "mean"),
  type = c("linear.prediction", "probability"),
  include.resid.var = TRUE,
  returnSims = FALSE,
  seed = NULL,
  .parallel = FALSE,
  .paropts = NULL,
  fix.intercept.variance = FALSE,
  ignore.fixed.terms = NULL
) {
  #--- Argument handling -------------------------------------------------------
  if (missing(newdata)) {
    newdata <- merMod@frame
  }
  if (any(c("data.frame") != class(newdata))) {
    if (any(c("tbl_df", "tbl") %in% class(newdata))) {
      newdata <- as.data.frame(newdata)
      warning("newdata is tbl_df or tbl object from dplyr package and has been
              coerced to a data.frame")
    } else {
      newdata <- as.data.frame(newdata)
    }
  }

  predict.type <- match.arg(type, c("linear.prediction", "probability"),
                            several.ok = FALSE)
  stat.type    <- match.arg(stat, c("median", "mean"), several.ok = FALSE)
  which.eff    <- match.arg(which, c("full", "fixed", "random", "all"),
                            several.ok = FALSE)

  # Set random seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  } else if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }

  # Check for probability type with linear models
  merMod.devcomp <- getME(merMod, "devcomp")
  if (merMod.devcomp$dims[["GLMM"]] == 0 &&
      merMod.devcomp$dims[["NLMM"]] == 0 &&
      predict.type == "probability") {
    predict.type <- "linear.prediction"
    warning("Asking for predictions on the probability scale makes no sense, resetting predict.type to linear.prediction",
            call. = FALSE)
  }

  # Detect if an intercept is present in random effects
  # This is checked for every call to maintain backward compatibility with
  # the message about intercept detection
  groupNames <- names(ngrps(merMod))
  hasREIntercept <- FALSE
  for (j in groupNames) {
    if (!is.na(names(attr(VarCorr(merMod)[[j]], "stddev")["(Intercept)"]))) {
      hasREIntercept <- TRUE
      break
    }
  }
  if (!hasREIntercept) {
    if (fix.intercept.variance) {
      fix.intercept.variance <- FALSE
    }
    message("No intercept detected, setting fix.intercept.variance to FALSE")
  }

  # Check for intercept when fix.intercept.variance is still TRUE
  if (fix.intercept.variance) {
    # Check if intercept is in fixed effects
    if (!"(Intercept)" %in% names(fixef(merMod))) {
      warning("No fixed-effect intercept detected. Variance adjustment may be unreliable.")
    }
  }

  # Check GLMM support
  check_glmm_support <- function(merMod, include.resid.var) {
    devcomp <- getME(merMod, 'devcomp')
    
    if (devcomp$dims[["GLMM"]] == TRUE) {
      family <- merMod@resp$family$family
      supported_families <- c("binomial", "poisson", "Gamma")
      
      if (!(family %in% supported_families)) {
        msg <- paste0(
          "Prediction for ", family, " GLMMs is not supported.\n",
          "Supported families: binomial, poisson, Gamma\n",
          "Try include.resid.var = FALSE as a fallback."
        )
        stop(msg, call. = FALSE)
      }
      
      if (include.resid.var && family == "binomial") {
        warning(paste0(
          "For binomial GLMMs, include.resid.var = TRUE simulates from the\n",
          "conditional binomial distribution (n-trial binomial simulation).\n",
          "This is the theoretically correct approach.\n",
          "To get predictions without residual variance, set include.resid.var = FALSE."
        ), call. = FALSE)
      }
    }
  }
  
  check_glmm_support(merMod, include.resid.var)

  #--- Simulations ------------------------------------------------------------
  # Order matters for reproducibility: sigma -> random effects -> fixed effects
  sigma_vec <- simulate_residual_variance(merMod, n.sims)
  random_list <- simulate_random_effects(
    merMod,
    newdata,
    n.sims,
    .parallel = .parallel,
    seed = seed
  )
  fixed_mat <- simulate_fixed_effects(
    merMod,
    newdata,
    n.sims,
    ignore.fixed.terms = ignore.fixed.terms,
    fix.intercept.variance = fix.intercept.variance,
    .parallel = .parallel,
    seed = seed
  )

  # Extract family/link info for GLMM support
  devcomp <- getME(merMod, 'devcomp')
  family_info <- NULL
  weights <- NULL
  
  if (devcomp$dims[["GLMM"]] == TRUE) {
    family_info <- list(
      family = merMod@resp$family$family,
      link = merMod@resp$family$link
    )
    
    # Handle weights for binomial GLMMs
    if (family_info$family == "binomial") {
      # First check for explicit weights
      if (!is.null(merMod@frame$weights)) {
        weights <- as.numeric(merMod@frame$weights)
      } else if (any(grepl("cbind", names(merMod@frame)))) {
        # For cbind syntax, extract weights from response columns
        cbind_cols <- names(merMod@frame)[grepl("cbind", names(merMod@frame))]
        if (length(cbind_cols) == 1) {
          weights <- rowSums(merMod@frame[, cbind_cols, drop = FALSE])
        }
      }
    }
  }

  #--- Combine components ------------------------------------------------------
combined_result <- combine_components(
     fixed_mat = fixed_mat,
     random_list = random_list,
     sigma_vec = sigma_vec,
     include.resid.var = include.resid.var,
     which = which.eff,
     family = if (!is.null(family_info)) family_info$family else NULL,
     link = if (!is.null(family_info)) family_info$link else NULL,
     weights = weights,
     use.probability = (include.resid.var && !is.null(family_info) &&
                        family_info$family %in% c("binomial", "poisson", "Gamma") &&
                        predict.type == "probability")
   )

  # Handle which = "all" case separately
  if (which.eff == "all") {
    yhat_arr <- combined_result$yhat
    pi.comps <- combined_result$components
  } else {
    yhat_arr <- combined_result
    pi.comps <- NULL
  }

  #--- Summarise --------------------------------------------------------------
  # For GLMMs with include.resid.var=TRUE, simulations are in response scale
  # so skip linkinv. For LMMs or GLMMs without residual variance, use original logic.
  devcomp <- getME(merMod, "devcomp")
  is.glmm.with.response <- (devcomp$dims[["GLMM"]] == TRUE && include.resid.var)
  
  outs <- summarise_predictions(
    yhat_arr = yhat_arr,
    level = level,
    stat.type = stat.type,
    predict.type = predict.type,
    N = nrow(newdata),
    merMod = merMod,
    which.eff = which.eff,
    pi.comps = pi.comps,
    is.glmm.with.response = is.glmm.with.response
  )

  # Attach simulation results if requested
  if (returnSims) {
    if (which.eff == "all") {
      attr(outs, "sim.results") <- pi.comps
    } else {
      attr(outs, "sim.results") <- yhat_arr
    }
  }

  outs
}
