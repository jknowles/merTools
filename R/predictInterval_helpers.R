#' Simulate residual variance for merMod objects
#'
#' Internal helper that reproduces the sigma-drawing logic from the original
#' `predictInterval()` implementation. It returns a numeric vector of length
#' `n.sims` containing the simulated residual standard deviations.
#'
#' @param merMod A merMod object from lme4.
#' @param n.sims Number of simulation draws.
#' @return Numeric vector of length `n.sims`.
#' @keywords internal
simulate_residual_variance <- function(merMod, n.sims) {
  devcomp <- getME(merMod, 'devcomp')
  # Linear model (no GLMM/NLMM)
  if (devcomp$dims[["GLMM"]] == 0 && devcomp$dims[["NLMM"]] == 0) {
    sigma <- sqrt(
      1 / rgamma(
        n.sims,
        shape = 0.5 * residDF.merMod(merMod),
        rate = 0.5 * devcomp$cmp[["pwrss"]]
      )
    )
    return(sigma)
  }

  # Binomial GLMM with logit/probit link – residual variance fixed at 1
  if (devcomp$dims[["GLMM"]] == TRUE &&
      merMod@resp$family$family == "binomial" &&
      merMod@resp$family$link %in% c('logit', 'probit')) {
    return(rep(1, n.sims))
  }

  # Fallback for other GLMM/NLMM cases (not extensively tested)
  warning(
    "Prediction for NLMMs or GLMMs that are not mixed binomial regressions is not tested. Sigma set at 1.",
    call. = FALSE
  )
  rep(1, n.sims)
}

#' Simulate fixed-effect predictions
#'
#' Draws `n.sims` samples from the multivariate normal distribution defined by
#' the fixed-effect estimates and their variance-covariance matrix, then
#' multiplies by the model matrix to produce predictions. It also respects
#' `ignore.fixed.terms` and optional intercept variance adjustment.
#'
#' @param merMod A merMod object.
#' @param newdata Data frame containing the prediction covariates.
#' @param n.sims Number of simulations.
#' @param ignore.fixed.terms Vector of fixed-effect names or indices to treat as known.
#' @param fix.intercept.variance Logical flag for intercept variance correction.
#' @param .parallel Logical flag to enable parallel execution via foreach.
#' @return Matrix of dimensions `nrow(newdata)` x `n.sims` containing simulated
#'   fixed-effect predictions.
#' @keywords internal
simulate_fixed_effects <- function(
  merMod,
  newdata,
  n.sims,
  ignore.fixed.terms = NULL,
  fix.intercept.variance = FALSE,
  .parallel = FALSE
) {
  fe.tmp <- fixef(merMod)
  vcov.tmp <- as.matrix(vcov(merMod))

  # Adjust intercept variance if requested (mirrors original logic)
  if (fix.intercept.variance) {
    rr <- ranef(merMod, condVar = TRUE)
    intercept.var <- vcov.tmp[1, 1]
    groupsizes <- ngrps(merMod)
    for (j in names(groupsizes)) {
      reMatrix <- attr(rr[[j]], which = 'postVar')
      groupVar <- (attr(VarCorr(merMod)[[j]], 'stddev')['(Intercept)'])^2
      extraPrec <- 0
      for (eff in seq_len(dim(reMatrix)[3])) {
        term <- 1 / (reMatrix[1, 1, eff] + groupVar)
        if (term > 0) {
          extraPrec <- extraPrec + term
        } else {
          warning("fix.intercept.variance got negative precision; better turn it off.")
        }
      }
      intercept.var <- intercept.var - 1 / extraPrec
    }
    if (intercept.var < 0) {
      warning("fix.intercept.variance got negative variance; better turn it off.")
    }
    ratio <- intercept.var / vcov.tmp[1, 1]
    prec.tmp <- solve(vcov.tmp)
    prec.tmp[1, 1] <- prec.tmp[1, 1] / ratio
    vcov.tmp[1, ] <- vcov.tmp[1, ] * ratio
    vcov.tmp <- solve(prec.tmp, tol = 1e-50)
  }

  # Inflate precision for ignored fixed terms (effectively zero variance)
  if (!is.null(ignore.fixed.terms)) {
    prec.tmp <- solve(vcov.tmp)
    for (term in ignore.fixed.terms) {
      prec.tmp[term, term] <- prec.tmp[term, term] * 1e15
    }
    vcov.tmp <- solve(prec.tmp, tol = 1e-50)
  }

  # Simulate coefficients - use parallel for large n.sims
  if (n.sims > 2000 || .parallel) {
    if (.parallel) {
      setup_parallel()
    }
    i <- 1:n.sims
    fe_call <- as.call(c(list(quote(foreach::foreach), i = i, .combine = 'rbind')))
    fe <- eval(fe_call)
    betaSim <- foreach::`%dopar%`(
      fe,
      mvtnorm::rmvnorm(n = 1, mean = fe.tmp, sigma = vcov.tmp, method = 'chol')
    )
  } else {
    betaSim <- abind::abind(
      lapply(1:n.sims, function(x) {
        mvtnorm::rmvnorm(n = 1, mean = fe.tmp, sigma = vcov.tmp, method = 'chol')
      }),
      along = 1
    )
  }
  colnames(betaSim) <- names(fe.tmp)
  rownames(betaSim) <- 1:n.sims

  # Build model matrix for newdata
  newdata.modelMatrix <- buildModelMatrix(merMod, newdata = newdata, which = "fixed")

  # Pad betaSim if model matrix has more columns (random slopes not in fixed effects)
  if (ncol(newdata.modelMatrix) > ncol(betaSim)) {
    pad <- matrix(
      rep(0),
      nrow = nrow(betaSim),
      ncol = ncol(newdata.modelMatrix) - ncol(betaSim)
    )
    if (ncol(pad) > 0) {
      message("Fixed effect matrix has been padded with 0 coefficients
            for random slopes not included in the fixed effects and interaction terms.")
    }
    colnames(pad) <- setdiff(colnames(newdata.modelMatrix), colnames(betaSim))
    betaSim <- cbind(betaSim, pad)
  }

  # Align columns between model matrix and coefficient matrix
  keep <- intersect(colnames(newdata.modelMatrix), colnames(betaSim))
  newdata.modelMatrix <- newdata.modelMatrix[, keep, drop = FALSE]
  betaSim <- betaSim[, keep, drop = FALSE]

  # Compute predictions: X %*% t(beta) -> (N x n.sims)
  fixed.xb <- newdata.modelMatrix %*% t(betaSim)
  fixed.xb
}

#' Compute random effect predictions for a single grouping factor
#'
#' Internal helper that multiplies predictor data by simulated random effect
#' coefficients for each grouping level. Handles new levels not present in the
#' original model data by returning zeros for their random effect contributions.
#'
#' @param data data.frame with predictor columns and a grouping factor column.
#'   The last column must be the grouping factor.
#' @param coefs 3D array of simulated coefficients with dimensions
#'   (n.sims x terms x levels).
#' @param group character name of the grouping factor column in data.
#' @return matrix of predicted values with dimensions (nrow(data) x n.sims).
#' @keywords internal
tmp.pred <- function(data, coefs, group) {
  # Identify new levels not in the model data
  new.levels <- unique(
    as.character(data[, group])[
      !as.character(data[, group]) %in% dimnames(coefs)[[3]]
    ]
  )


  if (length(new.levels) > 0) {
    msg <- paste0(
      "     The following levels of ", group, " from newdata \n -- ",
      paste0(new.levels, collapse = ", "),
      " -- are not in the model data. \n",
      "     Currently, predictions for these values are based only on the \n",
      " fixed coefficients and the observation-level error."
    )
    warning(msg, call. = FALSE)
  }

  # Initialize output array

  yhatTmp <- array(data = NA, dim = c(nrow(data), dim(coefs)[1]))
  colIdx <- ncol(data) - 1  # Number of predictor columns (excluding group)

  # Handle coefficient dimension mismatch

  # When predictor columns exceed coefficient terms, expand coefficients
  colLL <- length(1:colIdx)
  if (colLL > dim(coefs)[2]) {
    coefs_new <- array(
      NA,
      dim = c(dim(coefs)[1], colLL, dim(coefs)[3])
    )
    dimnames(coefs_new)[c(1, 3)] <- dimnames(coefs)[c(1, 3)]
    dimnames(coefs_new)[[2]] <- rep(dimnames(coefs)[[2]], dim(coefs_new)[2])
    for (k in 1:colLL) {
      coefs_new[, k, 1:dim(coefs)[3]] <- coefs[, 1, 1:dim(coefs)[3]]
    }
    coefs <- coefs_new
  }

  # Compute predictions for each observation
  for (i in 1:nrow(data)) {
    lvl <- as.character(data[, group][i])
    if (!lvl %in% new.levels) {
      # Known level: multiply predictors by coefficients
      yhatTmp[i, ] <- as.numeric(data[i, 1:colIdx]) %*% t(coefs[, 1:colIdx, lvl])
    } else {
      # New level: zero out the random effect contribution
      yhatTmp[i, ] <- rep(0, colIdx) %*% t(coefs[, 1:colIdx, 1])
    }
  }

  rownames(yhatTmp) <- rownames(data)
  return(yhatTmp)
}

#' Simulate random‑effect contributions for all grouping factors
#'
#' For each random effect term the function draws `n.sims` samples from the
#' conditional multivariate normal distribution (BLUPs + post‑var). The result is a
#' named list where each element is an array of dimensions `nrow(newdata)` ×
#' `n.sims`.
#' @param merMod A merMod object.
#' @param newdata Data frame containing the prediction covariates.
#' @param n.sims Number of simulations.
#' @param .parallel Logical flag to enable parallel execution via foreach.
#' @return List of matrices (rows = observations, cols = simulations).
simulate_random_effects <- function(
  merMod,
  newdata,
  n.sims,
  .parallel = FALSE
) {
  rr <- ranef(merMod, condVar = TRUE)
  re.xb <- vector(getME(merMod, 'n_rfacs'), mode = 'list')
  names(re.xb) <- names(ngrps(merMod))
  newdata.modelMatrix <- buildModelMatrix(model = merMod, newdata = newdata)

  for (j in names(re.xb)) {
    reMeans <- as.matrix(rr[[j]])
    reMatrix <- attr(rr[[j]], which = 'postVar')
    # Determine levels to keep based on newdata
    if (j %in% names(newdata)) {
      obslvl <- unique(as.character(newdata[, j]))
      alllvl <- rownames(reMeans)
      keep <- intersect(obslvl, alllvl)
    } else {
      obslvl <- colnames(newdata.modelMatrix)
      alllvl <- rownames(reMeans)
      keep <- intersect(obslvl, alllvl)
    }
    if (length(keep) > 0 && !identical(keep, alllvl)) {
      reMeans <- reMeans[keep, , drop = FALSE]
      dimnames(reMatrix)[[3]] <- alllvl
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else if (length(keep) > 0 && identical(keep, alllvl)) {
      dimnames(reMatrix)[[3]] <- alllvl
      reMatrix <- reMatrix[, , keep, drop = FALSE]
    } else {
      reMeans <- reMeans[1, , drop = FALSE]
      reMatrix <- reMatrix[, , 1, drop = FALSE]
    }

    # Simulate draws for each level
    tmpList <- vector('list', nrow(reMeans))
    for (k in seq_len(nrow(reMeans))) {
      meanTmp <- reMeans[k, ]
      matrixTmp <- as.matrix(reMatrix[, , k])
      tmpList[[k]] <- mvtnorm::rmvnorm(
        n = n.sims,
        mean = meanTmp,
        sigma = matrixTmp,
        method = 'chol'
      )
    }
    REcoefs <- sapply(tmpList, identity, simplify = 'array')
    dimnames(REcoefs) <- list(1:n.sims,
                               attr(reMeans, 'dimnames')[[2]],
                               attr(reMeans, 'dimnames')[[1]])

    # Build temporary data frame for prediction as in original code
    if (j %in% names(newdata)) {
      newdata.modelMatrix <- as.matrix(newdata.modelMatrix)
      tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
      tmp <- tmp[, !duplicated(colnames(tmp))]
      keepCols <- names(tmp)[names(tmp) %in% colnames(REcoefs)]
      if (length(keepCols) == 0) {
        keepCols <- grep(dimnames(REcoefs)[[2]], names(tmp), value = TRUE)
      }
      if (length(keepCols) == 0) {
        tmp <- cbind(model.frame(subbars(formula(merMod)), data = newdata),
                     var = newdata[, j])
        keepCols <- grep(dimnames(REcoefs)[[2]], names(tmp), value = TRUE)
      }
      if (length(keepCols) == 0) {
        # Add in an intercept for RE purposes
        tmp <- cbind(as.data.frame(newdata.modelMatrix), var = newdata[, j])
        tmp <- tmp[, !duplicated(colnames(tmp))]
        tmp <- cbind(data.frame(1), tmp)
        names(tmp)[1] <- "(Intercept)"
        keepCols <- "(Intercept)"
      }
      tmp <- tmp[, c(keepCols, "var"), drop = FALSE]
      tmp[, "var"] <- as.character(tmp[, "var"])
      colnames(tmp)[which(names(tmp) == "var")] <- names(newdata[, j, drop = FALSE])
      if (all(grepl(":", keepCols))) {
        # Strip out the interaction after
        keepCols <- unique(gsub("(.*):.*", "\\1", keepCols))
      }
    } else {
      # If newdata.modelMatrix is still sparse at this point, we need to convert it safely
      if (is(newdata.modelMatrix, 'dgCMatrix')) {
        newdata.modelMatrix <- as.matrix(newdata.modelMatrix)
        tmp <- as.data.frame(newdata.modelMatrix)
      } else {
        tmp <- as.data.frame(newdata.modelMatrix)
      }
      tmp <- tmp[, !duplicated(colnames(tmp))]
      # Assign grouping variable based on which model matrix column is active
      tmp$var <- names(tmp[keep])[max.col(tmp[keep])]
      keepCols <- names(tmp)[names(tmp) %in% dimnames(REcoefs)[[2]]]
      tmp <- tmp[, c(keepCols, "var"), drop = FALSE]
      tmp[, "var"] <- as.character(tmp[, "var"])
      colnames(tmp)[which(names(tmp) == "var")] <- j
    }

    # Compute contribution matrix using the internal prediction helper
    # Use parallel execution for large datasets
    if (nrow(tmp) > 1000 || .parallel) {
      if (requireNamespace("foreach", quietly = TRUE)) {
        if (.parallel) {
          setup_parallel()
        }
        # Split data into chunks of 500 rows
        tmp2 <- split(tmp, (seq_len(nrow(tmp)) - 1) %/% 500)
        tmp2 <- tmp2[vapply(tmp2, function(x) nrow(x) > 0, logical(1))]

        fe_call <- as.call(c(
          list(quote(foreach::foreach), i = seq_along(tmp2), .combine = 'rbind')
        ))
        fe <- eval(fe_call)

        re.xb[[j]] <- foreach::`%dopar%`(
          fe,
          tmp.pred(data = tmp2[[i]], coefs = REcoefs[, keepCols, , drop = FALSE], group = j)
        )
      } else {
        warning("foreach package is unavailable, parallel computing not available")
        re.xb[[j]] <- tmp.pred(data = tmp, coefs = REcoefs[, keepCols, , drop = FALSE], group = j)
      }
    } else {
      re.xb[[j]] <- tmp.pred(data = tmp, coefs = REcoefs[, keepCols, , drop = FALSE], group = j)
    }
  }
  re.xb
}

#' Combine fixed, random and residual components into a prediction array
#'
#' @param fixed_mat Matrix of dimensions `nrow(newdata)` × `n.sims`.
#' @param random_list List returned by `simulate_random_effects()`.
#' @param sigma_vec Residual standard deviations (output of `simulate_residual_variance`).
#' @param include.resid.var Logical, whether to add residual variance.
#' @param which Character: "full", "fixed", "random" or "all".
#' @return Array `nrow(newdata)` × `n.sims` (or list when `which="all"`).
#' @keywords internal
combine_components <- function(
  fixed_mat,
  random_list,
  sigma_vec,
  include.resid.var = TRUE,
  which = c('full', 'fixed', 'random', 'all')
) {
  which <- match.arg(which)
  re.xb <- random_list
  re.xb$fixed <- fixed_mat

  if (which == 'fixed') {
    yhat <- re.xb[['fixed']]
    if (include.resid.var) {
      N <- nrow(yhat)
      yhat <- abind::abind(
        lapply(seq_len(length(sigma_vec)), function(x) rnorm(N, yhat[, x], sigma_vec[x])),
        along = 2
      )
    }
    return(yhat)
  }

  if (which == 'random') {
    re.xb[['fixed']] <- NULL
    yhat <- Reduce('+', re.xb)
    if (include.resid.var) {
      N <- nrow(yhat)
      yhat <- abind::abind(
        lapply(seq_len(length(sigma_vec)), function(x) rnorm(N, yhat[, x], sigma_vec[x])),
        along = 2
      )
    }
    return(yhat)
  }

  if (which == 'full') {
    yhat <- Reduce('+', re.xb)
    if (include.resid.var) {
      N <- nrow(yhat)
      yhat <- abind::abind(
        lapply(seq_len(length(sigma_vec)), function(x) rnorm(N, yhat[, x], sigma_vec[x])),
        along = 2
      )
    }
    return(yhat)
  }

  if (which == 'all') {
    yhat <- Reduce('+', re.xb)
    # For 'all', apply residual variance to components (for returnSims)
    # and also to yhat for the summary
    if (include.resid.var) {
      N <- nrow(yhat)
      # Apply residual variance to each component
      for (i in seq_along(re.xb)) {
        re.xb[[i]] <- abind::abind(
          lapply(seq_len(length(sigma_vec)), function(x) rnorm(N, re.xb[[i]][, x], sigma_vec[x])),
          along = 2
        )
      }
      # Apply residual variance to combined yhat
      yhat <- abind::abind(
        lapply(seq_len(length(sigma_vec)), function(x) rnorm(N, yhat[, x], sigma_vec[x])),
        along = 2
      )
    }
    return(list(yhat = yhat, components = re.xb))
  }
}

#' Summarise prediction simulations into intervals
#'
#' Takes the simulated prediction array and computes summary statistics
#' (fit, lower bound, upper bound) based on the specified statistic type
#' and confidence level. Handles probability transformation and the
#' `which = "all"` case with component-wise summaries.
#'
#' @param yhat_arr matrix of simulated predictions with dimensions (N x n.sims),
#'   or for `which = "all"`, the combined yhat from `combine_components()`.
#' @param level numeric confidence level (e.g., 0.8 for 80
#' @param stat.type character either "median" or "mean" for the central estimate.
#' @param predict.type character either "linear.prediction" or "probability".
#' @param N integer number of observations in the original newdata.
#' @param merMod merMod object, required when `predict.type = "probability"` to
#'   access the link function.
#' @param which.eff character the effect type: "full", "fixed", "random", or "all".
#' @param pi.comps list of component prediction matrices (from `combine_components()`
#'   when `which = "all"`), or NULL for other cases.
#' @return data.frame with columns `fit`, `upr`, `lwr`. For `which = "all"`,
#'   includes additional columns `effect` and `obs`.
#' @keywords internal
summarise_predictions <- function(
  yhat_arr,
  level,
  stat.type,
  predict.type,
  N,
  merMod = NULL,
  which.eff = "full",
  pi.comps = NULL
) {
  # Calculate quantile bounds

upCI <- 1 - ((1 - level) / 2)
  loCI <- (1 - level) / 2

  # Initialize output data frame
  outs <- data.frame(
    fit = rep(NA_real_, N),
    upr = rep(NA_real_, N),
    lwr = rep(NA_real_, N)
  )

  # Compute summary statistics based on stat.type
  if (stat.type == "median") {
    outs[, 1:3] <- t(apply(yhat_arr, 1, quantile,
                           prob = c(0.5, upCI, loCI),
                           na.rm = TRUE))
  } else if (stat.type == "mean") {
    outs$fit <- apply(yhat_arr, 1, mean, na.rm = TRUE)
    outs[, 2:3] <- t(apply(yhat_arr, 1, quantile,
                           prob = c(upCI, loCI),
                           na.rm = TRUE))
  }

  # Apply link function for probability predictions
  if (predict.type == "probability") {
    if (is.null(merMod)) {
      stop("merMod object required for probability predictions", call. = FALSE)
    }
    if (nrow(outs) == 1) {
      outs <- t(apply(outs, 2, merMod@resp$family$linkinv))
    } else {
      outs <- apply(outs, 2, merMod@resp$family$linkinv)
    }
    outs <- as.data.frame(outs)
    names(outs) <- c("fit", "upr", "lwr")
  }

  # Handle which = "all" case with component breakdown
  if (which.eff == "all" && !is.null(pi.comps)) {
    # Process each component
    for (i in seq_along(pi.comps)) {
      if (stat.type == "median") {
        pi.comps[[i]] <- t(apply(pi.comps[[i]], 1, quantile,
                                 prob = c(0.5, upCI, loCI),
                                 na.rm = TRUE))
        pi.comps[[i]] <- as.data.frame(pi.comps[[i]])
        names(pi.comps[[i]]) <- c("fit", "upr", "lwr")
      } else if (stat.type == "mean") {
        tmp <- pi.comps[[i]]
        pi.comps[[i]] <- data.frame(
          fit = rep(NA_real_, N),
          upr = NA_real_,
          lwr = NA_real_
        )
        pi.comps[[i]]$fit <- apply(tmp, 1, mean, na.rm = TRUE)
        pi.comps[[i]][, 2:3] <- t(apply(tmp, 1, quantile,
                                        prob = c(upCI, loCI),
                                        na.rm = TRUE))
      }
      # Apply link function for probability predictions
      if (predict.type == "probability") {
        pi.comps[[i]] <- apply(pi.comps[[i]], 2, merMod@resp$family$linkinv)
        pi.comps[[i]] <- as.data.frame(pi.comps[[i]])
        names(pi.comps[[i]]) <- c("fit", "upr", "lwr")
      }
    }

    # Bind components together with effect labels
    componentOut <- dplyr::bind_rows(pi.comps, .id = "effect")
    outs <- cbind(data.frame(effect = "combined"), outs)
    outs <- suppressWarnings(dplyr::bind_rows(outs, componentOut))
    outs$obs <- rep(seq_len(N), nrow(outs) %/% N)
  }

  as.data.frame(outs)
}
