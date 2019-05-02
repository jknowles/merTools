#' Display model fit summary of x or x like objects, fast
#'
#' Faster than the implementation in the arm package because it avoids refitting
#'
#' @title fastdisp: faster display of model summaries
#' @param x a model object
#' @param ... additional arguments to pass to \code{arm::\link[arm]{display}}
#' including number of digits
#' @details The time saving is only noticeable for large, time-consuming (g)lmer
#' fits.
#' @import arm
#' @return A printed summary of a x object
#' @examples
#' \donttest{
#' #Compare the time for displaying this modest model
#' require(arm)
#' m1 <- lmer(y ~ lectage + studage + (1|d) + (1|s), data=InstEval)
#' system.time(display(m1))
#' system.time(fastdisp(m1))
#' }
#' @seealso \code{\link[arm]{display}}
#' @rdname fastdisp
#' @export fastdisp
fastdisp <- function (x, ...) {
  UseMethod("fastdisp", x)
}

#' @rdname fastdisp
#' @importFrom stats df pt
#' @export
fastdisp.merMod <- function (x, ...)
{
  .local <- function (x, digits = 2, detail = FALSE)
  {
    out <- NULL
    out$call <- x@call
    print(out$call)
    fcoef <- fixef(x)
    useScale <- getME(x, "devcomp")$dims["useSc"]
    corF <- vcov(x)@factors$correlation
    coefs <- cbind(fcoef, corF@sd)
    if (length(fcoef) > 0) {
      if (!useScale) {
        coefs <- coefs[, 1:2, drop = FALSE]
        out$z.value <- coefs[, 1]/coefs[, 2]
        out$p.value <- 2 * pnorm(abs(out$z.value), lower.tail = FALSE)
        coefs <- cbind(coefs, `z value` = out$z.value,
                       `Pr(>|z|)` = out$p.value)
      }
      else {
        out$t.value <- coefs[, 1]/coefs[, 2]
        coefs <- cbind(coefs, `t value` = out$t.value)
      }
      dimnames(coefs)[[2]][1:2] <- c("coef.est", "coef.se")
      if (detail) {
        pfround(coefs, digits)
      }
      else {
        pfround(coefs[, 1:2], digits)
      }
    }
    out$coef <- coefs[, "coef.est"]
    out$se <- coefs[, "coef.se"]
    cat("\nError terms:\n")
    vc <- easyVarCorr(VarCorr(x), useScale = useScale,
                            digits)
    print(vc[, c(1:2, 4:ncol(vc))], quote = FALSE)
    out$ngrps <- lapply(x@flist, function(x) length(levels(x)))
    is_REML <- isREML(x)
    llik <- logLik(x, REML = is_REML)
    out$AIC <- AIC(llik)
    # out$deviance <- deviance(refitML(x))
    out$n <- getME(x, "devcomp")$dims["n"]
    # Dhat <- -2 * (llik)
    # pD <- out$deviance - Dhat
    # out$DIC <- out$deviance + pD
    cat("---\n")
    cat(sprintf("number of obs: %d, groups: ", out$n))
    cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "),
              collapse = "; "))
    cat(sprintf("\nAIC = %g", round(out$AIC, 1)))
    # cat(round(out$DIC, 1))
    # cat("\ndeviance =", fround(out$deviance, 1), "\n")
    if (useScale < 0) {
      out$sigma.hat <- sigma(x)
      cat("overdispersion parameter =", fround(out$sigma.hat,
                                               1), "\n")
    }
    return(invisible(out))
  }
  .local(x, ...)
}

#' @rdname fastdisp
#' @export
fastdisp.merModList <- function(x, ...){
  .local <- function (x, digits = 2, detail = FALSE)
  {
  out <- NULL
  useScale <- getME(x[[1]], "devcomp")$dims["useSc"]
  #useScale <- TRUE
  out$call <- x[[1]]@call
  listFE <- modelFixedEff(x)
  row.names(listFE) <- listFE$term
  out$t.value <- listFE$statistic
  out$coef <- listFE$estimate
  out$se <- listFE$std.error
  listRE <- modelRandEffStats(x)
  out$ngrps <- lapply(x[[1]]@flist, function(x) length(levels(x)))
  is_REML <- isREML(x[[1]])
  llik <- lapply(x, logLik, REML = is_REML)
  out$AIC <- mean(unlist(lapply(llik, AIC)))
  out$n <- round(mean(unlist(lapply(lapply(lapply(x, getME, "devcomp"),
                                     "[[", "dims"), "[", 2))), 0) # round to nearest integer
  print(out$call)
  if (!detail) {
    pfround(listFE[, 2:3], digits)
  }
  else {
    listFE$p.value <- 2 * pt(abs(listFE$statistic), listFE$df, lower.tail = FALSE)
    pfround(listFE[, 2:6], digits)
  }
  cat("\nError terms:\n")
  vc <- easyVarCorr(VarCorr(x[[1]]), useScale = useScale,
                    digits)
  # Resort the output of the random effect summary
  listRE <- listRE[grep("cor_", listRE$term, invert=TRUE), ]
  resid <- listRE[listRE$group == "Residual", ]
  listRE <- listRE[listRE$group != "Residual", ]
  listRE <- rbind(listRE, resid)
  #
  vc[, 3] <- as.character(round(listRE$estimate^2, digits = digits))
  vc[, 4] <- as.character(round(listRE$estimate, digits = digits))
  print(vc[, c(1:2, 4:ncol(vc))], quote = FALSE)
  cat("---\n")
  cat(sprintf("number of obs: %d, groups: ", out$n))
  cat(paste(paste(names(out$ngrps), out$ngrps, sep = ", "),
            collapse = "; "))
  cat(sprintf("\nAIC = %g", round(out$AIC, 1)))
  cat("---\n")
  # cat(round(out$DIC, 1))
  # cat("\ndeviance =", fround(out$deviance, 1), "\n")
  if (useScale < 0) {
    out$sigma.hat <- sigma(x)
    cat("overdispersion parameter =", fround(out$sigma.hat,
                                             1), "\n")
    cat("---\n")
  }
  return(invisible(out))
  }
  .local(x, ...)
}

